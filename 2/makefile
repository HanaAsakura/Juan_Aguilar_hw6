burger.gif:burger_data.dat burger.py 
	python burger.py
burger_data.dat: burger.x
	./burger.x > burger_data.dat
burger.x: burger.c
	cc burger.c -o burger.x -lm
