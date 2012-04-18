cluster: *.cpp *.h
	g++ *.cpp *.h -o cluster -O3 -march=core2
clean:
	rm cluster
