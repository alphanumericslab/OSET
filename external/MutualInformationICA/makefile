CFLAGS=g++
OPTIONS= -O3


milca: miutils.C miutils.h milcadelay.C ICAtests.C MIhigherdim.C MIxnyn.C MIClustering.C
	$(CFLAGS) -c miutils.C -o miutils.o $(OPTIONS)
	$(CFLAGS) milca.C -o milca miutils.o $(OPTIONS)
	$(CFLAGS) milcadelay.C -o milcadelay miutils.o $(OPTIONS)
	$(CFLAGS) ICAtests.C -o ICAtests miutils.o $(OPTIONS)
	$(CFLAGS) MIhigherdim.C -o MIhigherdim miutils.o $(OPTIONS)
	$(CFLAGS) MIxnyn.C -o MIxnyn miutils.o $(OPTIONS)
	$(CFLAGS) MIClustering.C -o MIClustering miutils.o $(OPTIONS)

