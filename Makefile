cluster: *.cpp *.h maze
	g++ *.cpp *.h -o cluster
maze: clique/*.h clique/*.c
	make -C clique
clean:
	rm cluster clique/mace
