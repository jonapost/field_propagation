name := testNTST
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

CPPFLAGS += -g
CPPFLAGS += -O3
CPPFLAGS += -msse
CPPFLAGS += -ftree-vectorizer-verbose=3
CPPFLAGS += -ffast-math


ifdef LHEP
  CPPFLAGS += -DflagLHEP
endif

ifdef QGSP
  CPPFLAGS += -DflagQGSP
endif

ifdef QGSP_EMV
  CPPFLAGS += -DflagQGSP_EMV
endif

ifdef QGSC
  CPPFLAGS += -DflagQGSC
endif

ifdef FTFP
  CPPFLAGS += -DflagFTFP
endif

ifdef QGSP_BIC
  CPPFLAGS += -DflagQGSP_BIC
endif

ifdef QGSP_BERT
  CPPFLAGS += -DflagQGSP_BERT
endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk



