CFLAGS=g++
OPTIONS= -O3


milca: miutils.c miutils.h milcadelay.c ICAtests.c MIhigherdim.c MIxnyn.c MIClustering.c
	$(CFLAGS) -c miutils.c -o miutils.o $(OPTIONS)
	$(CFLAGS) milca.c -o milca miutils.o $(OPTIONS)
	$(CFLAGS) milcadelay.c -o milcadelay miutils.o $(OPTIONS)
	$(CFLAGS) ICAtests.c -o ICAtests miutils.o $(OPTIONS)
	$(CFLAGS) MIhigherdim.c -o MIhigherdim miutils.o $(OPTIONS)
	$(CFLAGS) MIxnyn.c -o MIxnyn miutils.o $(OPTIONS)
	$(CFLAGS) MIClustering.c -o MIClustering miutils.o $(OPTIONS)

