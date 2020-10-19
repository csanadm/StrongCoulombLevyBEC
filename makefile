CXX   = g++
LD    = g++
SRCSUF  = cpp
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit -lMinuit2

WrkDir     := $(shell pwd)
ExeDir     = exe
ObjDir     = object
DepDir     = object/deps
SrcDir     = sources
ClsDir     = classes

CLASSES  = Levy_reader HypSInt_reader CorrFunc_reader CorrFunc_evaluator StrongCorrFunc_evaluator StrongCorrFunc_reader
PROGRAMS = AfterBurner_CorrFunc_bettereta Levy_filler HypSInt_filler CorrFunc_filler 
CPPS     = functions plotter $(PROGRAMS)

CLASS_SOURCES = $(addprefix $(ClsDir)/,$(addsuffix .$(SRCSUF), $(CLASSES)))
PROGRAM_SOURCES = $(addprefix $(SrcDir)/,$(addsuffix .$(SRCSUF), $(CPPS)))

SOURCES =
SOURCES += $(PROGRAM_SOURCES)
SOURCES += $(CLASS_SOURCES)
ALL_SOURCES = $(sort $(SOURCES))

define COMPILE_TEMPLATE
-include $(DepDir)/$(notdir $(basename $(1))).d
$(ObjDir)/$(notdir $(basename $(1))).o:
	$(CXX) -c $(ROOTCFLAGS) -I$(WrkDir)/$(SrcDir) -I$(WrkDir)/$(ClsDir) -MD -MP -MF $(DepDir)/$(notdir $(basename $(1))).d $(1) -o $$@ 
endef
$(foreach source, $(ALL_SOURCES), $(eval $(call COMPILE_TEMPLATE,$(source))))

exe/Levy_filler.exe: object/Levy_filler.o object/functions.o object/plotter.o
	$(LD) -O -m32 $^ $(ROOTLIBS) -o $@

exe/HypSInt_filler.exe: object/HypSInt_filler.o object/functions.o
	$(LD) -O -m32 $^ $(ROOTLIBS) -o $@

exe/AfterBurner_CorrFunc_bettereta.exe: object/AfterBurner_CorrFunc_bettereta.o object/functions.o
	$(LD) -O -m32 $^ $(ROOTLIBS) -o $@

exe/CorrFunc_filler.exe: object/CorrFunc_filler.o object/HypSInt_reader.o object/Levy_reader.o object/functions.o
	$(LD) -O -m32 $^ $(ROOTLIBS) -o $@

all:   	$(addprefix exe/,$(addsuffix .exe, $(PROGRAMS)))

clean:
	rm -f $(ExeDir)/*.* $(ObjDir)/*.* $(DepDir)/*.*
