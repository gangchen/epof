cluster: *.cpp *.h mace
	g++ *.cpp *.h -o cluster
mace: clique/*.h clique/*.c
	make -C clique
clean:
	rm cluster clique/mace
