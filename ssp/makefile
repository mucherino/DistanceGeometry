
# SSP project

OBJ= main.o ssp.o bp.o utils.o

ssp: $(OBJ)
	g++ -O3 -o ssp $(OBJ) -lm

.c.o:
	gcc -O3 -c $<

clean:
	\rm -f *.o ssp

