cluster: *.cpp *.h mace
	g++ *.cpp *.h -o cluster -pg
mace: clique/*.h clique/*.c
	make -C clique
clean:
	rm cluster clique/mace
