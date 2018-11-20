CXXFLAGS=-Wall 
PROGS=EF2d-main-Perio_static EF2d-main-Perio_inst EF2d-main-Perio_move
OBJS=CG.o		EF2d-base.o
LIBS = -lm
all: $(PROGS)
	echo fait! 
EF2d-main-Perio_static:EF2d-main-Perio_static.o $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLGAS) $(LDFLAGS) $(LIBS) 
EF2d-main-Perio_inst:EF2d-main-Perio_inst.o $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLGAS) $(LDFLAGS) $(LIBS)
EF2d-main-Perio_move:EF2d-main-Perio_move.o $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLGAS) $(LDFLAGS) $(LIBS)

clean:
	-rm *.o $(PROGS) *~ 
EF2d-main-Perio_static.o: CG.hpp EF2d-base.hpp R2.hpp clock.hpp
EF2d-main-Perio_inst.o: CG.hpp EF2d-base.hpp R2.hpp clock.hpp 
EF2d-main-Perio_move.o: CG.hpp EF2d-base.hpp R2.hpp clock.hpp 
EF2d-base.o:EF2d-base.hpp R2.hpp clock.hpp
