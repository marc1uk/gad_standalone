#all: main monitor draw
all: comp

main: checkgad.cpp
	g++ -std=c++11 -g $< -o $@ `root-config --cflags --libs` -lMinuit

monitor: monitor_gad_calibration.cpp
	g++ -g -Wno-psabi $< -o $@ `root-config --cflags --libs`

draw: just_draw.cpp
	g++ -g $< -o $@ `root-config --cflags --libs`

comp: compare_shapes.cpp
	g++ -g $< -o $@ `root-config --cflags --libs`
