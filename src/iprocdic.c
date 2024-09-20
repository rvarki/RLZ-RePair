#include <stdlib.h>
#include <stdio.h>


// iprocdic.c
//
// Intermediate tool for bigRepair for integers 
// derived from procdic.c and modified for the case
// in which the input alphabet of prefix free parsing consists of integers

// transform the PFP dictionary adding a 
// a unique terminator after each string
// the values used for unique terminators are U, U+1, and so on ...
// where U is currently 0x78000000 (constant Unique in makefile)
// the program check that input symbols are indeed smaller than U 
// otherwise it exits with an error message

// The constant Unique is defined in the makefile to make
// sure the same value is used also in ipostproc.
// Unique is the first symbol used as a separator; all input symbols must 
// be smaller than it. Unique can be made larger, but all unique 
// separators must be smaller than 2^{31}  (this is checked in the code)
 
// This version uses the new (Spire'19) dictionary format where 
// the lengths (in symbols) of dictionary strings are provided in a 
// separate .len file of int32_t integers. 


int main (int argc, char **argv)
{ 
  char foname[1024], flname[1024];
  FILE *fi,*fo, *fl;
  int n=Unique;  // first integer value used as a separator
  int c;
  fputs("==== Command line:\n",stderr);
  for(int i=0;i<argc;i++)
    fprintf(stderr," %s",argv[i]);
  fputs("\n",stderr);

  if(argc!=2) {
     fprintf(stderr,"Usage: %s <dictionaryfilename>\n\n", argv[0]);
     exit(1);
  }

  // open dictionary file
  fi = fopen(argv[1],"r");
  if (fi == NULL) { 
    fprintf(stderr,"Cannot open file %s\n",argv[1]);
    exit(1);
  }
  // open dictionary lengths file
  sprintf(flname,"%s.len",argv[1]);
  fl = fopen(flname,"r");
  if (fl == NULL) { 
    fprintf(stderr,"Cannot open file %s\n",flname);
    exit(1);
  }
  // open integer output file
  sprintf(foname,"%s.int",argv[1]);
  fo = fopen(foname,"w");
  if (fo == NULL) { 
    fprintf(stderr,"Cannot create file %s\n",foname);
    exit(1);
  }
  
  // main loop 
  while (1) {
    // read length of next word 
    int wlen;
    int e = fread(&wlen,sizeof(int),1,fl);
    if(e!=1) {
      if(feof(fl)) break; // end file
      else {perror("Error reading dictionary word length"); exit(1); }
    }
    // read actual word
    for(int i=0;i<wlen;i++) {
      e = fread(&c,sizeof(int),1,fi);
      if(e!=1) {perror("Unexpected end of dictionary"); exit(1); }
      if(c<0 || c>=Unique) {fprintf(stderr, "Dictionary symbol %x larger than %x\n",c,Unique-1); exit(1);} 
      e = fwrite (&c,sizeof(int),1,fo); // write int to output file 
      if(e!=1) {perror("Error writing to .int file"); exit(1); }
    }
    e = fwrite (&n,sizeof(int),1,fo); // write terminator as an int to output file 
    if(e!=1) {perror("Error writing to .int file"); exit(1); }
    if(++n < 0) {// update separator and check it is still valid 
      fprintf(stderr,"No more valid separator symbol. Exiting...\n"); exit(1);
    }
  }
  c = getc(fi); // there should be no further chars
  if(c!=EOF) {perror("Unexpected trailing chars in dictionary"); exit(1); }
  fclose(fo);
  fclose(fi);
  fprintf (stderr,"%i strings\n",n-Unique);
  fputs("=== Preprocessing completed!\n",stderr);
  return 0;
}
