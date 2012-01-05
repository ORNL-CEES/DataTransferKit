MOAB_BASE_DIR = /home/stuart/software/builds/moab

include ${MOAB_BASE_DIR}/lib/iMesh-Defs.inc

#CXXFLAGS = -O3 -march=native -Wall
CXXFLAGS = -g -Wall		
LDFLAGS = ${IMESH_LIBS} -L${MOAB_BASE_DIR}/lib
CXX = g++
CXXSOURCES = coupler/Transfer_Evaluator.cc coupler/Transfer_Map.cc coupler/Data_Transfer_Manager.cc example/Post_Proc.cc example/Physics_A.cc example/Physics_B.cc example/TE_Physics_A.cc example/TE_Physics_B.cc example/driver.cc 
CXXOBJS = ${CXXSOURCES:.cc=.o}
EXECUTABLE = simple_transfer

${EXECUTABLE}: ${CXXOBJS} Makefile
	${CXX} ${CXXFLAGS} -o $@ ${CXXOBJS} ${LDFLAGS}

Transfer_Evaluator.o: coupler/Transfer_Evaluator.cc coupler/Transfer_Evaluator.hh

Transfer_Map.o: coupler/Transfer_Map.cc coupler/Transfer_Map.hh coupler/Transfer_Evaluator.hh

Data_Transfer_Manager.o: coupler/Data_Transfer_Manager.cc coupler/Data_Transfer_Manager.hh coupler/Transfer_Evaluator.hh coupler/Transfer_Map.hh

Post_Proc.o: example/Post_Proc.cc example/Post_Proc.hh

Physics_A.o: example/Physics_A.cc example/Physics_A.hh example/Post_Proc.hh

Physics_B.o: example/Physics_B.cc example/Physics_B.hh example/Post_Proc.hh

TE_Physics_A.o: example/TE_Physics_A.cc example/TE_Physics_A.hh coupler/Transfer_Evaluator.hh

TE_Physics_B.o: example/TE_Physics_B.cc example/TE_Physics_B.hh coupler/Transfer_Evaluator.hh

driver.o: example/driver.cc example/Physics_A.hh example/Physics_B.hh example/TE_Physics_A.hh example/TE_Physics_B.hh coupler/Transfer_Evaluator.hh coupler/Data_Transfer_Manager.hh

.cc.o:
	${CXX} ${CXXFLAGS} ${IMESH_CXXFLAGS} ${IMESH_INCLUDES} -o $@ -c $<

clean:
	 rm -f example/*.o coupler/*.o ${EXECUTABLE}