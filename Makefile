all: main.c Lab4_IO.c
	mpicc -o main main.c Lab4_IO.c -lm

serial:
	gcc -o serialtester serialtester.c Lab4_IO.c -lm

trim:
	gcc -o trim datatrim.c Lab4_IO.c -lm

clean:
	rm datagen main serialtester
