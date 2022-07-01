all: main monitor draw

main: checkgad.cpp
	g++ -std=c++11 -g $< -o $@ `root-config --cflags --libs` -lMinuit

monitor: monitor_gad_calibration.cpp
	g++ -g $< -o $@ `root-config --cflags --libs`

draw: just_draw.cpp
	g++ -g $< -o $@ `root-config --cflags --libs`
