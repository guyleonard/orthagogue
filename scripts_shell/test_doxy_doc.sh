#!/bin/sh
# Create *.png files from text @startuml code and output *.png images to ./doc-src/images folder
# Recusively search from current folder scanning *.c, *.cpp, *.h, and *.doc files
java  -Djava.awt.headless=true -jar $PLANTUML_JAR -v -o $PWD/doc-src/images  "./**.(c|cpp|doc|h)"
doxygen