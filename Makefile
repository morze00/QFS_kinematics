CXX= clang++
LINKER= clang++
LDFLAGS= `root-config --libs`
CPPFLAGS= `root-config --cflags`
CXXFLAGS= -Wno-deprecated -g
LNKFLAGS= -g
EXEC=qfs

all: $(EXEC)

OBJ=  qfs.o

HEADERS=  headers.hh	\
	info.hh	\
	tree.hh	\
	libs.hh

%.o : %.cc $(DEPS)
	@$(MAKEDEPEND)
	@${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
	@echo "	CXX $@"

%.o : %.C $(DEPS)
	@$(MAKEDEPEND)
	@${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
	@echo "	CXX $@"

$(EXEC): $(OBJ)
	@${LINKER} ${LDFLAGS} ${LNKFLAGS} -o $(EXEC) $(OBJ)
	@echo "	COMP $(EXEC)"

clean:
	rm -f $(OBJ)
	rm -f $(EXEC)
