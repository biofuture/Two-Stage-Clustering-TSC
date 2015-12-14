OBJECTS = global.o  kmerneedle.o herarchicalcluster.o  greedheuristic.o tsc_main.o
TSC:$(OBJECTS)
	g++ -o TSC_version_1.2_i86linux64 -O3  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -fstack-check -fvisibility=hidden -pipe -fomit-frame-pointer -mthreads -m64 -g -fopenmp $(OBJECTS)
global.o: common.h global.cpp global.h
	g++ -O 3 -c global.cpp 
kmerneedle.o: kmerneedle.h kmerneedle.cpp
	g++ -I /home/jiang/bin/seqan -fopenmp -O3 -c kmerneedle.cpp
herarchicalcluster.o: herarchicalcluster.h herarchicalcluster.cpp
	g++ -O3 -c herarchicalcluster.cpp
greedheuristic.o: greedheuristic.h greedheuristic.cpp
	g++ -O3 -c greedheuristic.cpp
tsc_main.o: global.h common.h tsc_main.cpp
	g++ -O3 -c  tsc_main.cpp
.PHONY : clean
clean :
	rm  (TSC $(objects))
