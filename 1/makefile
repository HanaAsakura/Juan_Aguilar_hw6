E:=1
ALPHA:=1

particle_in_field.pdf: trayectoria_$(E)_$(ALPHA).dat
	python particle_in_field.py trayectoria_$(E)_$(ALPHA).dat
trayectoria_$(E)_$(ALPHA).dat: particle_in_field.x
	./particle_in_field.x $(E) $(ALPHA)
particle_in_field.x: particle_in_field.c
	cc particle_in_field.c -o particle_in_field.x -lm

