COMPILLER=g++
FLAGS=-Wextra -Wall -Werror -pedantic -std=c++11

all: start

start: main.o
	$(COMPILLER) $(FLAGS) -o solution main.o

main.o: main.cpp
	$(COMPILLER) -c $(FLAGS) main.cpp

linear: linear.o
	$(COMPILLER) $(FLAGS) -o linear linear.o

linear.o: linear.cpp
	$(COMPILLER) -c $(FLAGS) linear.cpp

clean:
	@-rm -f *.o *.gch *.dat solution linear
	@echo "Clean success"
