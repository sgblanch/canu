CC = gcc
CXX = g++
CFLAGS := -O4 -fexpensive-optimizations -finline-functions -fomit-frame-pointer -pipe 
CXXFLAGS := $(CFLAGS)
LOBJS = ssw.o
LCPPOBJS = ssw_cpp.o
LIB = libssw.so
EXAMPLE_CPP = example_cpp

.PHONY: all default clean

default: $(EXAMPLE_CPP) $(LIB) 

all: default

java: $(JAVA_JAR) $(JAVA_LIB)

$(LIB): ssw.C ssw.H
	$(CC) $(CFLAGS) -fPIC -shared  -o $@ $<

$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) -lm -lz

$(EXAMPLE_CPP): example.C $(LOBJS) $(LCPPOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) -lm -lz

$(JAVA_LIB): sswjni.c ssw.c ssw.h
	$(CC) $(CFLAGS) $(JAVA_INLCUDES) -fPIC -shared -o $@ $< ssw.c 

$(JAVA_JAR): $(JAVA_OBJ)
	jar cvfe $@ ssw.Example $^

%.class: %.java
	javac $<
	
ssw.o: ssw.C ssw.H
	$(CC) -c -o $@ $< $(CFLAGS)

ssw_cpp.o: ssw_cpp.C ssw_cpp.H ssw.H
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	-rm -f $(LOBJS) $(LCPPOBJS) $(PROG) $(LIB) $(EXAMPLE) $(EXAMPLE_CPP) $(JAVA_LIB) $(JAVA_JAR) $(JAVA_OBJ) *~ 
