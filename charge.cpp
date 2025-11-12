#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void charge(int me, int n, int np, int *iBeg, int *iEnd) {
    int r= n%np;

    if (me<r)
    {
    *iBeg = me * (n / np)+me;
    *iEnd = (me + 1) * (n / np) + me;
    }
    else{
    *iBeg = r*(n/np+1)+(me-r)*(n/np);
    *iEnd = *iBeg+(n/np)-1;
    }
    
}
