CC=icpc
CCFLAGS = -c -O3 -fopenmp -xMIC-AVX512 -ltbb 
LDFLAGS = -fopenmp  -ltbb 
SOURCES = main.cc data_manager.cc filter.cc
OBJECTS = $(SOURCES:.cc=.o)
EXECUTABLE = PUNAS

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

main.o: main.cc
	$(CC) $(CCFLAGS) main.cc -o main.o

data_manager.o: data_manager.cc data_manager.h
	$(CC) $(CCFLAGS) data_manager.cc -o data_manager.o

filter.o: filter.cc filter.h
	$(CC) $(CCFLAGS) filter.cc -o filter.o
	
clean:
		rm -f $(EXECUTABLE) *.o
