#all: main monitor draw comp
all: main

#main: checkgad.cpp
main: gad_standalone.cpp
	g++ -std=c++11 -g -Wno-psabi $< -o $@ `root-config --cflags --libs` -lMinuit

monitor: monitor_gad_calibration.cpp
	g++ -g -Wno-psabi $< -o $@ `root-config --cflags --libs`

draw: just_draw.cpp
	g++ -g -Wno-psabi $< -o $@ `root-config --cflags --libs`

comp: compare_shapes.cpp
	g++ -g -Wno-psabi $< -o $@ `root-config --cflags --libs`
