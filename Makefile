cluster: *.cpp *.h mace
	g++ *.cpp *.h -o cluster -O3 -march=core2
mace: clique/*.h clique/*.c
	make -C clique
clean:
	rm cluster clique/mace
