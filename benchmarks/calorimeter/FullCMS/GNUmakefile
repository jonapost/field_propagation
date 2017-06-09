name := full_cms
G4TARGET := $(name)
G4EXLIB := true

include $(G4INSTALL)/config/architecture.gmk
 
.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifndef XERCESCROOT	
	@echo XERCESCROOT not defined!
endif
