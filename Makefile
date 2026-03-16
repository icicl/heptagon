OBJS = objs/heptagon.o objs/node.o objs/heptagon_cu.o

# link .o -> bin
heptagon: $(OBJS)
	nvcc $(OBJS) -o heptagon -lgomp

# compile .c -> .o
objs/heptagon.o: heptagon.c
	mkdir -p objs
	gcc -fopenmp -O3 -g -c $< -o $@
objs/node.o: node.c
	mkdir -p objs
	gcc -O3 -g -c $< -o $@

# compile .cu -> .o
objs/heptagon_cu.o: heptagon.cu
	mkdir -p objs
	nvcc -O3 -g -c $< -o $@
