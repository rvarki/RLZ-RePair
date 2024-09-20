#include <stdlib.h>
#include <stdio.h>

// transform a dictionary converting chars to int32_t's and adding a 
// a unique terminator after each string
// the values used for unique terminators are 256, 257, and so on ... 
// this version uses the new (Spire'19) dictionary format where 
// the lengths of dictionary strings are provided in a separate .len file
// in int32_t integers. 

// input: file.dicz & file.dicz.len
// output: file.dicz.int

int main (int argc, char **argv)
{ 
  char foname[1024], flname[1024];
  FILE *fi,*fo, *fl;
  int n=256;  // first integer value used as a separator
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
      c = getc(fi); // next char
      if(c==EOF) {perror("Unexpected end of dictionary"); exit(1); }
      e = fwrite (&c,sizeof(int),1,fo); // write char as an int to output file 
      if(e!=1) {perror("Error writing to .int file"); exit(1); }
    }
    e = fwrite (&n,sizeof(int),1,fo); // write terminator as an int to output file 
    if(e!=1) {perror("Error writing to .int file"); exit(1); }
    n++; // update terminator 
  }
  c = getc(fi); // there should be no further chars
  if(c!=EOF) {perror("Unexpected trailing chars in dictionary"); exit(1); }
  fclose(fo);
  fclose(fi);
  fprintf(stderr,"%i strings\n",n-256);
  fputs("=== Preprocessing completed!\n",stderr);
  exit(0);
}
