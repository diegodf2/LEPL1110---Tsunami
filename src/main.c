# include "tsunami.h"

int main(void)
{   

        char *meshName = "../data/PacificMedium.txt";
        char *resultBaseName = "../output/tsunamiMedium";

        double dt = 4;
        int nMax  = 4000;
        int sub   = 100;

        tsunamiCompute(dt,nMax,sub,meshName,resultBaseName);
        tsunamiAnimate(dt,nMax,sub,meshName,resultBaseName);
         

        exit(0);     
}