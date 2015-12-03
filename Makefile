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

exe: lp mkl_test lp-mkl


lp-mkl: lp-mkl.o mt19937p.o
	$(CC) $(CFLAGS) $(OPTFLAGS) $(LDFLAGS) $^ -o $@ $(LIBMKL) $(INCMKL)

lp-mkl.o: lp-mkl.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(LDFLAGS) $< $(LIBMKL) $(INCMKL)


lp: lp.o mt19937p.o
	$(CC) $(CFLAGS) $(OPTFLAGS) $(LDFLAGS) $^ -o $@

lp.o: lp.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(LDFLAGS) $<

mkl_test: mkl_test.o 
	$(CC)  $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBMKL) $(INCMKL)

mkl_test.o: mkl_test.c
	$(CC) -c $(CFLAGS) $(LDFLAGS)  $< $(LIBMKL) $(INCMKL)

%.o: %.c
	$(CC) -c $(CFLAGS) $(LDFLAGS) $<

# === Cleanup and tarball

clean:
	rm -f *.o

realclean: clean
	rm -f lp
