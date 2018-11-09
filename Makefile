CXX = g++
CXXFLAGS = -g -O0

SRCDIR = src
BUILDDIR = build
TARGET = qmc

SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS = $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(BUILDDIR)/%.o)
INC = -I include

all: qmc

$(TARGET): $(OBJECTS)
	@echo "Linking...."
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)";
	$(RM) -r $(BUILDDIR) $(TARGET)

cleandocs:
	-rm -fr html/ latex/

.PHONY : clean
