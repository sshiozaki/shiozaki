COMPILER = icpc

#local environment
DEST	 = ~/bin
LIBDIR = /home/shiozaki/lib
TARDIR = ./bin
OBJDIR = ./obj
SRCDIR = ./src

CPPFLAGS = -Wall -fast -std=c++0x -fopenmp -MMD -MP -Wextra -Winit-self -Wno-unused-parameter -Wfloat-equal
LDFLAGS  = -fopenmp -lifcore
LIBS     = 
INCLUDE  = -I./include

#TextParser
TXTP      = $(LIBDIR)/TextParser-1.8.5
LDFLAGS += -L$(TXTP)/lib -lTP
INCLUDE += -I$(TXTP)/include

#SDFlib
SDFLib = $(LIBDIR)/SDF
LDFLAGS += -L$(SDFLib)/lib -lSDF
INCLUDE += -I$(SDFLib)/include/sdf

#hdf5
HDF5 = $(LIBDIR)/hdf5-1.8.20
LDFLAGS += -L$(HDF5)/lib -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -ldl -Wl,-rpath -Wl,$(HDF5)/lib
INCLUDE += -I$(HDF5)/include

TARGET  = $(TARDIR)/FDstentDeployment

ifeq "$(strip $(SRCDIR))" ""
  SRCDIR = .
endif
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR = .
endif

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))
DEPENDS = $(OBJECTS:.o=d)

$(TARGET): $(OBJECTS) $(LIBS)
	@[ -d $(TARDIR) ] || mkdir -p $(TARDIR)
	$(COMPILER) -o $@ $^ $(LDFLAGS) 

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(COMPILER) $(CPPFLAGS) $(INCLUDE) -o $@ -c $<


all: clean $(TARGET)

clean:
	rm -f $(OBJECTS) $(DEPENDS) $(TARGET)

install:$(TARGET)
	install -s $(TARGET) $(DEST)

-include $(DEPENDS)

