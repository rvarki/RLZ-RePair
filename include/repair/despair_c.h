#ifndef DESPAIR_H
#define DESPAIR_H

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

int despair(const char *base_name);

#ifdef __cplusplus
}
#endif

#endif