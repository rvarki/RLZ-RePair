EXECS = irepair idespair repair despair procdic postproc iprocdic ipostproc
CFLAGS = -Wall -g -m64 -O9

# these excutables are not used by BigRepair 
#   idespair repair despair toint 
# they compile but have not been tested!

all: $(EXECS)

# targets not producing a file declared phony
.PHONY: all clean


procdic: procdic.c
	gcc $(CFLAGS) -o $@ $< -D Unique=0x78000000

postproc: postproc.c
	gcc $(CFLAGS) -o $@ $< -D Unique=0x78000000

iprocdic: iprocdic.c
	gcc $(CFLAGS) -o $@ $< -D Unique=0x78000000

ipostproc: ipostproc.c
	gcc $(CFLAGS) -o $@ $< -D Unique=0x78000000

repair: repair.o array.o hash.o heap.o records.o basics.o makefile
	gcc $(CFLAGS) -o repair repair.o array.o hash.o heap.o records.o basics.o

repair.o: repair.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c repair.c

despair: despair.o basics.o makefile
	gcc $(CFLAGS) -o despair despair.o basics.o

despair.o: despair.c basics.h makefile
	gcc $(CFLAGS) -c despair.c

irepair: irepair.o array.o hash.o heap.o records.o basics.o makefile
	gcc $(CFLAGS) -o irepair irepair.o array.o hash.o heap.o records.o basics.o

irepair.o: irepair.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c irepair.c

idespair: idespair.o basics.o makefile
	gcc $(CFLAGS) -o idespair idespair.o basics.o

idespair.o: idespair.c basics.h makefile
	gcc $(CFLAGS) -c idespair.c

array.o: array.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c array.c

hash.o: hash.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c hash.c

heap.o: heap.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c heap.c

records.o: records.c array.h hash.h heap.h records.h basics.h makefile
	gcc $(CFLAGS) -c records.c

basics.o: basics.c basics.h makefile
	gcc $(CFLAGS) -c basics.c

toint: toint.c makefile
	gcc $(CFLAGS) -o toint toint.c

clean:
	rm -f $(EXECS) *.o

