CXX = mpic++
CXXFLAGS = -Wall -std=c++11

all: parsplit

parsplit: parsplit.c
	$(CXX) $(CXXFLAGS) --prefix /usr/local/share/OpenMPI -o $@ $<

clean:
	rm -f parsplit