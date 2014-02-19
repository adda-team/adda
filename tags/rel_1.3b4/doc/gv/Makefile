# Processes GraphWiz files in the same directory using the dot tool
SHELL := /bin/sh
GVSRC := $(wildcard *.gv)
EXT   := svg

# Tests if dot is available directly, otherwise uses a hard-coded path 
ifeq ($(shell which dot > /dev/null 2>&1 && echo 0),0)
  DOT := dot
else
  DOT := "C:\Program Files (x86)\Graphviz2.30\bin\dot"
endif
# Relative path to folder with wiki files (images are copied there on install) 
WIKIPATH := ../../web/wiki

OUTPUT   := $(GVSRC:.gv=.$(EXT))

.PHONY: all clean

all: $(OUTPUT)

$(OUTPUT): %.$(EXT): %.gv
	$(DOT) -T$(EXT) -o$@ $<

clean:
	rm *.$(EXT) 

install: all
	cp *.$(EXT) $(WIKIPATH)
