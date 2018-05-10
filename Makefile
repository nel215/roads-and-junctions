all:
	g++ -DLOCAL --std=c++0x -W -Wall -Wno-sign-compare -O2 -s -pipe -mmmx -msse -msse2 -msse3 main.cpp
RoadsAndJunctionsVis.class: ./RoadsAndJunctionsVis.java
	javac ./RoadsAndJunctionsVis.java
