#ifndef REPAIR_C
#define REPAIR_C

#ifdef __cplusplus
extern "C" {
#endif

#include "basics.h"
#include "hash.h"
#include "heap.h"
#include "records.h"
#include "bitarr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

int repair_main(const char *seqfile, const char *rlzFile, const char *refFile, float factor, int MB);

#ifdef __cplusplus
}
#endif

#endif