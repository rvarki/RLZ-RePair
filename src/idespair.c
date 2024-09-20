// decompress an int-based .R .C pair

/*

Repair -- an implementation of Larsson and Moffat's compression and
decompression algorithms.
Copyright (C) 2010-current_year Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of
Chile. Blanco Encalada 2120, Santiago, Chile. gnavarro@dcc.uchile.cl

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include "basics.h"

relong u; // |text| and later current |C| with gaps

int *C; // compressed text

int c;  // real |C|

int alph; // max used terminal symbol

Tpair *R; // rules

int n; // |R|

char *ff;
FILE *f;
int maxdepth = 0;

// expand symbol i at depth d
// if i is a terminal it is writen to f
// if is a non terminal recusively expand left and right children 
int expand (int i, int d)

   { int ret = 1;
     while (i >= alph) // while i is a non-terminal
       { ret += expand(R[i-alph].left,d+1);
         i = R[i-alph].right; d++; // second recursive call replaced by iteration
       }
     if (fwrite(&i,sizeof(int),1,f) != 1)
  { fprintf (stderr,"Error: cannot write file %s\n",ff);
    exit(1);
  }
     if (d > maxdepth) maxdepth = d;
     return ret;
   }

int main (int argc, char **argv)

   { char fname[PATH_MAX], outname[PATH_MAX];
     FILE *Tf,*Rf,*Cf;
     int i,len,c,u;
     struct stat s;
     fputs("==== Command line:\n",stderr);
     for(int i=0;i<argc;i++)
       fprintf(stderr," %s",argv[i]);
     fputs("\n",stderr);     
     if (argc != 2)
  { fprintf (stderr,"Usage: %s <filename>\n"
        "Decompresses <filename> from its .C and .R extensions.\n"
        "Decompressed file is <filename>.out\n"
        "This is a version for integer sequences\n",argv[0]);
    exit(1);
  }
     strcpy(fname,argv[1]);
     strcat(fname,".R");
     if (stat (fname,&s) != 0)
  { fprintf (stderr,"Error: cannot stat file %s\n",fname);
    exit(1);
  }
     len = s.st_size;
     Rf = fopen (fname,"r");
     if (Rf == NULL)
  { fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }
  // read alphabet size for the original input string 
     if (fread(&alph,sizeof(int),1,Rf) != 1)
  { fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
     n = (len-sizeof(int))/sizeof(Tpair); // number of rules 
     R = (void*)malloc(n*sizeof(Tpair));  // array of rules 
     if (fread(R,sizeof(Tpair),n,Rf) != n)
  { fprintf (stderr,"Error: cannot read file %s\n",fname);
    exit(1);
  }
     fclose(Rf);

     // open C file and get the number of symbols in it 
     strcpy(fname,argv[1]);
     strcat(fname,".C");
     if (stat (fname,&s) != 0)
  { fprintf (stderr,"Error: cannot stat file %s\n",fname);
    exit(1);
  }
     c = len = s.st_size/sizeof(int); // number of symbols in C string
     Cf = fopen (fname,"r");
     if (Cf == NULL)
  { fprintf (stderr,"Error: cannot open file %s for reading\n",fname);
    exit(1);
  }
  
     // open output file
     strcpy(outname,argv[1]);
     strcat(outname,".out");
     Tf = fopen (outname,"w");
     if (Tf == NULL)
  { fprintf (stderr,"Error: cannot open file %s for writing\n",outname);
    exit(1);
  }
     u = 0; f = Tf; ff = outname;
     // read one symbol of C at a time and expand it 
     for (;len>0;len--)
  { if (fread(&i,sizeof(int),1,Cf) != 1)
       { fprintf (stderr,"Error: cannot read file %s\n",fname);
         exit(1);
       }
    u += expand(i,0);
  }
     fclose(Cf);
     if (fclose(Tf) != 0)
  { fprintf (stderr,"Error: cannot close file %s\n",outname);
    exit(1);
  }
  // here n is the number of rules, n+alpha the effective alphabet in C 
  long est_size = (long) ( (2.0*n+(n+c)*(float)blog(n+alph-1))/8) + 1;  
  fprintf (stderr,"IDesPair succeeded\n");
  fprintf (stderr,"   Original ints: %i\n",u);
  fprintf (stderr,"   Size of the original input alphabet: %i\n",alph);
  fprintf (stderr,"   Number of rules: %i\n",n);
  fprintf (stderr,"   Compressed sequence length: %i\n",c);
  fprintf (stderr,"   Maximum rule depth: %i\n",maxdepth);
  fprintf (stderr,"   Estimated output size (bytes): %ld\n",est_size);
  fprintf (stderr,"   Compression ratio: %0.2f%%\n", (100.0*8*est_size)/(u*blog(alph-1)));
  // original estimate:   ((4.0*(n-alph)+((n-alph)+c)*(float)blog(n-1))/(u*blog(alph-1))*100.0));
  exit(0);
   }

