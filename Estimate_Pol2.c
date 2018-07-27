#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include <time.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>

#define row 10001
#define column 100001

///// matrix.c /////
#ifdef __cplusplus

template<typename T>T**AllocMatrix(int u,int v)
{
    int i; T**a,*b;
    try { a=(T**)new char[(sizeof*a+sizeof*b*v)*u]; }
    catch (...) { a=0; }
    if (a) b=(T*)(a+u); else return 0;
    for (i=0;i<u;i++,b+=v) a[i]=b;
    return a;
}
#define ALLOC_MATRIX(T,U,V) AllocMatrix<T>(U,V)
#define FREE(X) delete[]X

#else

void*AllocMatrix(int s,int u,int v)
{
    int i,t=s*v; char**a,*b;

    a=(char**)malloc((sizeof*a+t)*u);
    if (a) b=(char*)(a+u); else return 0;
    for (i=0;i<u;i++,b+=t) a[i]=b;
    return a;
}
#define ALLOC_MATRIX(T,U,V) (T**)AllocMatrix(sizeof(T),U,V)
#define FREE(X) free(X)

#endif
///// /////
