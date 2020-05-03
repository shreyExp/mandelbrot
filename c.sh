#!/bin/bash
g++ -g -o $1 main.cpp `pkg-config --libs --cflags opencv4 jsoncpp`
