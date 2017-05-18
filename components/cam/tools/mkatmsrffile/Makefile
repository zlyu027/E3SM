CASEROOT=/home/bhillma/projects/cmdv-rrm/ACME/cime/cases/skybridge/FC5AQUAP.enax4v1
include $(CASEROOT)/Macros.make

#LIBDIR=/ptmp/jedwards/SMS_D.f09_g16.B1850CN.bluefire.C.083815/lib
LIBDIR=/gscratch/bhillma/acme_scratch/skybridge/FC5AQUAP.enax4v1/bld/lib
LIBS = -L$(LIBDIR) -lcsm_share -lpio -lmct -lmpeu -L$(LIB_NETCDF) -lnetcdf -L$(LIB_PNETCDF) -lpnetcdf
MODDIR=$(LIBDIR)/include


mkatmsrffile: mkatmsrffile.F90 
	$(SFC) $(FFLAGS) -I$(MODDIR) $< -o $@  $(LIBS) $(LDFLAGS) 


