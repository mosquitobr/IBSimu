CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

tofl203d: tofl203d.o
	$(CC) -o tofl203d tofl203d.o $(LDFLAGS)
tofl203d.o: tofl203d.cpp

clean:
	$(RM) *.o tofl203d 

