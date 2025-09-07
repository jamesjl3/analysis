AUTOMAKE_OPTIONS = foreign subdir-objects

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -I/sphenix/user/jamesj3j3/local/RooUnfold/include \
  -I$(srcdir)/JetValidation \
  -I$(srcdir)/Jets_miniPAC \
  -I$(srcdir)/JetResoStudy \
  -I$(srcdir)/EMJetVal \
  -I$(srcdir)/EEC_PAC \
  -I$(srcdir)/JetCorrection \
  -I$(srcdir)/JetUnfolding \
  -I$(srcdir)/RooUnfoldJets

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L/sphenix/user/jamesj3j3/local/RooUnfold/lib

pkginclude_HEADERS = \
  JetValidation/JetValidation.h \
  Jets_miniPAC/Jets_miniPAC.h \
  JetResoStudy/JetResoStudy.h \
  EMJetVal/EMJetValUnfolded.h \
  JetResoStudy/JetResoStudy_SUBT.h \
  EEC_PAC/EEC_Reso.h \
  JetCorrection/JetCorrection.h \
  JetUnfolding/JetUnfolding.h \
  JetMatchingSubjets/JetMatchingSubjets.h \
  JetMatchingSubjets/ZsjTools.h \
  RooUnfoldJets/RooUnfoldJets.h \
  JetUnfoldingSubjets.h \
  JetUnfoldingSubjets_noUniquePtr.h

lib_LTLIBRARIES = \
  libJetValidation.la

ROOTDICTS = \
  JetUnfoldingSubjets_noUniquePtr_LinkDef.h


pcmdir = $(libdir)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libJetValidation_la_SOURCES = \
  JetValidation/JetValidation.cc \
  Jets_miniPAC/Jets_miniPAC.cc \
  JetResoStudy/JetResoStudy.cc \
  EMJetVal/EMJetValUnfolded.cc \
  JetResoStudy/JetResoStudy_SUBT.cc \
  EEC_PAC/EEC_Reso.cc \
  JetCorrection/JetCorrection.cc \
  JetUnfolding/JetUnfolding.cc \
  RooUnfoldJets/RooUnfoldJets.cc \
  JetMatchingSubjets/JetMatchingSubjets.cc \
  JetMatchingSubjets/ZsjTools.cc \	
  JetUnfoldingSubjets.cc \
  JetUnfoldingSubjets_noUniquePtr.cc

libJetValidation_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lg4dst \
  -lg4eval \
  -lg4jets \
  -ltrackbase_historic \
  -lRooUnfold

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libJetValidation.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
