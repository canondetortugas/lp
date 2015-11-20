#
# To build with a different compiler / on a different platform, use
#     make PLATFORM=xxx
#
# where xxx is
#     icc = Intel compilers
#     gcc = GNU compilers
#     clang = Clang compiler (OS X default)
#
# Or create a Makefile.in.xxx of your own!
#

PLATFORM=icc
include Makefile.in.$(PLATFORM)

.PHONY: exe clean realclean

# === Executables

exe: lp

lp: lp.o 
	$(CC) $(OMP_CFLAGS) $^ -o $@

path.o: lp.c
	$(CC) -c $(OMP_CFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

# === Cleanup and tarball

clean:
	rm -f *.o

realclean: clean
	rm -f lp
