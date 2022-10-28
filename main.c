#include "main.h"

int main()
{
    
Coeffs = makeLowShelf (48000, 100, 1, 10);
    
double fq[26] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,20000};
double out[26] = {0};
    
getMagnitudeForFrequencyArray(&fq, &out,26, 48000, Coeffs);
    
    for(int i =0; i<26; i++)
    {
     printf("FQ: %f, PH: %f \n", fq[i], out[i]);
    }
return (0);
}
