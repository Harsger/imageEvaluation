all:matchNfit

matchNfit: matchNfit.C
	g++ `pkg-config opencv --cflags` matchNfit.C -o matchNfit `pkg-config opencv --libs` `root-config --glibs --cflags` -lSpectrum -std=c++0x

clean:
	rm matchNfit
