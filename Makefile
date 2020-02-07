all:matchNfit clip

matchNfit: matchNfit.C
	g++ `pkg-config opencv --cflags` matchNfit.C -o matchNfit `pkg-config opencv --libs` `root-config --glibs --cflags` -lSpectrum -std=c++0x

clip: clip.C
	g++ `pkg-config opencv --cflags` clip.C -o clip `pkg-config opencv --libs` `root-config --glibs --cflags` -lSpectrum -std=c++0x

clean:
	rm matchNfit clip
