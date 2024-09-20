#include <iostream>

extern "C" float factor = 0.75f;

extern "C" {
    #include "heap.h"
    #include "array.h"
    #include "basics.h"
    #include "records.h"
    #include "hash.h"
}


int main()
{
    Trarray Rec = createRecords(factor,10);
    Tpair mypair;
    mypair.left = 10;
    mypair.right = 100;
    insertRecord(&Rec, mypair);
    Theap Heap = createHeap(100,&Rec,factor,20);
    std::cout << "Hello World" << std::endl; 
    return 0;
}