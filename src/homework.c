# include "tsunami.h"

typedef struct {
    double *X;
    double *Y;
    double *bath;
    int nNode;
    int nElem;
    int *elem;
    double *E;
    double *U;
    double *V;
    double *FE;
    double *FU;
    double *FV;
    double *USol;
    double *VSol;
    double *ESol;
    int nEdge;
    int *Edges;
} myStruct;

typedef struct {
    double *K1E;
    double *K2E;
    double *SumE;
    double *K1U;
    double *K2U;
    double *SumU;
    double *K1V;
    double *K2V;
    double *SumV;
    double dt;
} HeunStruct;



void freeEverything(myStruct *tsunami, HeunStruct *Heun) {
    free(Heun->K1U);
    free(Heun->K2U);
    free(Heun->SumU);
    free(Heun->K1E);
    free(Heun->K2E);
    free(Heun->SumE);
    free(Heun->K1V);
    free(Heun->K2V);
    free(Heun->SumV);
    free(Heun);
    free(tsunami->bath);
    free(tsunami->elem);
    free(tsunami->X);
    free(tsunami->Y);
    free(tsunami->Edges);
    free(tsunami->USol);
    free(tsunami->VSol);
    free(tsunami->ESol);
    free(tsunami);
}

void createStructureBasedOnFile(const char *meshFileName, myStruct *tsunami) {
    int i,j,nNode,nElem,nEdge,trash;
    
    FILE* file = fopen(meshFileName,"r");

    fscanf(file, "Number of nodes %d \n",&nNode); 
    double *X    = malloc(sizeof(double)*nNode);
    double *Y    = malloc(sizeof(double)*nNode);  
    double *bath = malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&X[i],&Y[i],&bath[i]); 

    fscanf(file, "Number of triangles %d \n",&nElem); 
    int *elem = malloc(sizeof(int)*3*nElem);
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);   

    int *edges = malloc(sizeof(int)*nEdge*4);
    fscanf(file,"Number of edges %d \n",&nEdge);
    for (i = 0; i < nEdge; i++) {
        fscanf(file,"%d : %d %d : %d %d \n", &trash, &edges[4*i], &edges[4*i+1], &edges[4*i+2], &edges[4*i+3]);
    }
    fclose(file); 


    tsunami->nEdge  = nEdge;
    tsunami->bath   = bath;
    tsunami->elem   = elem;
    tsunami->nElem  = nElem;
    tsunami->X      = X;
    tsunami->Y      = Y;
    tsunami->nNode  = nNode;
    tsunami->Edges  = edges;
    tsunami->USol   = malloc(sizeof(double) * (nElem*3+1));
    tsunami->VSol   = malloc(sizeof(double) * (nElem*3+1));
    for (i = 0; i < nElem*3+1; i++) {
        tsunami->USol[i] = 0.;
        tsunami->VSol[i] = 0.;
    }
    
}

void mapElem(myStruct *tsunami, int indice, int mapTriangle[3]) {
    for (int i = 0; i < 3; i++)
        mapTriangle[i] = indice*3 + i;
}

void computedphi(double x[3], double y[3], double phi_x[3], double phi_y[3], double J) {
    phi_x[0] = (y[1] - y[2])/J;
    phi_x[1] = (y[2] - y[0])/J;
    phi_x[2] = (y[0] - y[1])/J;
    phi_y[0] = (x[2] - x[1])/J;
    phi_y[1] = (x[0] - x[2])/J;
    phi_y[2] = (x[1] - x[0])/J;
}

void integrateOnTriangles(myStruct *tsunami) {
    double *X       = tsunami->X;
    double *Y       = tsunami->Y;
    int *elem       = tsunami->elem;
    int nElem       = tsunami->nElem;
    double *U       = tsunami->U;
    double *FU      = tsunami->FU;
    double *E       = tsunami->E;
    double *FE      = tsunami->FE;
    double *V       = tsunami->V;
    double *FV      = tsunami->FV;
    double *bath    = tsunami->bath;


    int i,j,k,mapTriangle[3];
    double phi_x[3],phi_y[3],xLoc[3],yLoc[3],u,v,e,h,f,J,x,y;
    double phi[3][3] = {{1-gaussTriangleXsi[0]-gaussTriangleEta[0], gaussTriangleXsi[0], gaussTriangleEta[0]},
                        {1-gaussTriangleXsi[1]-gaussTriangleEta[1], gaussTriangleXsi[1], gaussTriangleEta[1]},
                        {1-gaussTriangleXsi[2]-gaussTriangleEta[2], gaussTriangleXsi[2], gaussTriangleEta[2]}};
    
    for (i = 0; i < nElem; i++) {
        mapElem(tsunami, i, mapTriangle);
         
        for (j = 0; j < 3; j++) {
            xLoc[j] = X[elem[i*3 + j]];
            yLoc[j] = Y[elem[i*3 + j]];
        }
        J = fabs((xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]));
        computedphi(xLoc, yLoc, phi_x, phi_y, J);
        for(k = 0; k < 3; k++) {
            u = 0; e = 0; v = 0; h = 0; x = 0; y = 0;
            for (j = 0; j < 3; j++) {
                x += phi[k][j] * xLoc[j];
                y += phi[k][j] * yLoc[j];
                u += phi[k][j] * U[mapTriangle[j]];
                v += phi[k][j] * V[mapTriangle[j]];
                e += phi[k][j] * E[mapTriangle[j]];
                h += phi[k][j] * bath[elem[i*3+j]];
            }
            double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  	        double lat = asin(z3d/R)*180./PI;
            f = 2.*Omega * sin(lat);
            for (j = 0; j < 3; j++) {
                FE[mapTriangle[j]] += ((phi_x[j]*h*u + phi_y[j]*h*v) * (4*R*R + x*x + y*y) / (4*R*R) + phi[k][j] * (h*(x*u + y*v)) / (R*R)) * J * gaussTriangleWeight[k];
                FU[mapTriangle[j]] += (phi[k][j] * (f*v - Gamma*u) + phi_x[j]*g*e * (4*R*R + x*x + y*y) / (4*R*R) + phi[k][j] * g*x*e / (2*R*R)) * J * gaussTriangleWeight[k];
                FV[mapTriangle[j]] += (-phi[k][j] * (f*u+Gamma*v) + phi_y[j] *g*e * (4*R*R + x*x + y*y) / (4*R*R) + phi[k][j] * g*y*e / (2*R*R)) * J * gaussTriangleWeight[k];
            }
        }
    }
}

void mapEdge(myStruct *tsunami, int indice, int mapEdge[2][2]) {
    int *node         = &tsunami->Edges[4*indice];
    int *elem         = &tsunami->Edges[4*indice+2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mapEdge[j][i] = tsunami->nElem * 3;
            if (elem[j] != -1) {
                for (int k = 0; k < 3; k++) {
                    if (tsunami->elem[elem[j]*3 + k] == node[i]) {
                        mapEdge[j][i] = elem[j]*3 + k;
                    }
                }
            } 
        }
    }
}

void normal(double x_xsi, double y_xsi, double n[3]) {
    n[2] = sqrt(x_xsi*x_xsi + y_xsi*y_xsi);
    n[0] = y_xsi / n[2];
    n[1] = -x_xsi / n[2];
}

void integrateOnEdges(myStruct *tsunami) {
    double *X       = tsunami->X;
    double *Y       = tsunami->Y;
    int *elem       = tsunami->elem;
    int nElem       = tsunami->nElem;
    int nEdges      = tsunami->nEdge;
    double *U       = tsunami->U;
    double *FU      = tsunami->FU;
    double *E       = tsunami->E;
    double *FE      = tsunami->FE;
    double *V       = tsunami->V;
    double *FV      = tsunami->FV;
    double *bath    = tsunami->bath;
    int *edges      = tsunami->Edges;
    
    int i,j,k,map[2][2];
    double xLoc[2],yLoc[2],x_xsi,y_xsi,norm[3],J;
    double eG,eD,uG,uD,vG,vD,h,unG,unD,estar,ustar,x,y;

    double phi[2][2] = {{(1. - gaussEdgeXsi[0]) / 2., (1. + gaussEdgeXsi[0]) / 2.},
                        {(1. - gaussEdgeXsi[1]) / 2., (1. + gaussEdgeXsi[1]) / 2.}};
    
    for (i = 0; i < nEdges; i++) {
        mapEdge(tsunami, i, map);
        for (j = 0; j < 2; j++) {
            xLoc[j] = X[edges[4*i+j]];
            yLoc[j] = Y[edges[4*i+j]];
        }
        
        int border = (map[1][0] == 3*nElem);

        x_xsi = xLoc[1] - xLoc[0];
        y_xsi = yLoc[1] - yLoc[0];
        
        normal(x_xsi,y_xsi,norm);
        J = norm[2] / 2.;
        
        for (j = 0; j < 2; j++) {
            eG = 0.; eD = 0.; uG = 0.; uD = 0.; vG = 0.; vD = 0.;
            x = 0.; y = 0.; h = 0.;
            for (k = 0; k < 2; k++) {
                x  += phi[j][k] * xLoc[k];
                y  += phi[j][k] * yLoc[k];
                h  += phi[j][k] * bath[edges[4*i+k]];
                eG += phi[j][k] * E[map[0][k]];
                eD += phi[j][k] * E[map[1][k]];
                uG += phi[j][k] * U[map[0][k]];
                uD += phi[j][k] * U[map[1][k]];
                vG += phi[j][k] * V[map[0][k]];
                vD += phi[j][k] * V[map[1][k]];
            }
        
            if (border) eD = eG;
            unG = uG*norm[0] + vG*norm[1];
            unD = border ? -unG : uD*norm[0] + vD*norm[1];
            estar = (eG + eD) / 2. + sqrt(h / g) * (unG - unD) / 2.;
            ustar = (unG + unD) / 2. + sqrt(g / h) * (eG - eD) / 2.;

            for (k = 0; k < 2; k++) {
                FE[map[0][k]] -= phi[j][k] * h * ustar * (4*R*R + x*x + y*y) / (4*R*R) * J * gaussEdgeWeight[j];
                FE[map[1][k]] += phi[j][k] * h * ustar * (4*R*R + x*x + y*y) / (4*R*R) * J * gaussEdgeWeight[j];
                FU[map[0][k]] -= phi[j][k] * norm[0] * estar * g * (4*R*R + x*x + y*y) / (4*R*R) * J * gaussEdgeWeight[j];
                FU[map[1][k]] += phi[j][k] * norm[0] * estar * g * (4*R*R + x*x + y*y) / (4*R*R) * J * gaussEdgeWeight[j];
                FV[map[0][k]] -= phi[j][k] * norm[1] * estar * g * (4*R*R + x*x + y*y) / (4*R*R) * J * gaussEdgeWeight[j];
                FV[map[1][k]] += phi[j][k] * norm[1] * estar * g * (4*R*R + x*x + y*y) / (4*R*R) * J * gaussEdgeWeight[j];
            }
        }
    }
}

void Phi_iPhiJOnEdgeAndInverse(myStruct *tsunami) {
    double *FE  = tsunami->FE;
    double *FU  = tsunami->FU;
    double *FV  = tsunami->FV;
    double *X   = tsunami->X;
    double *Y   = tsunami->Y;
    int *elem   = tsunami->elem;
    int nElem   = tsunami->nElem;

    int i,j,k,map[3];

    double xLoc[3],yLoc[3],J,tempE[3], tempU[3], tempV[3];

    double inverse[3][3] = {{18.,-6.,-6.},
                            {-6.,18.,-6.},
                            {-6.,-6.,18.}};

    for (i = 0; i < nElem; i++) {
        mapElem(tsunami, i, map);
        for (j = 0; j < 3; j++) {
            xLoc[j] = X[elem[i*3 + j]];
            yLoc[j] = Y[elem[i*3 + j]];
        }
        J = fabs((xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]));

        for (j = 0; j < 3; j++) {
            tempE[j]    = FE[map[j]];
            FE[map[j]]  = 0;
            tempU[j]    = FU[map[j]];
            FU[map[j]]  = 0;
            tempV[j]    = FV[map[j]];
            FV[map[j]]  = 0;
        }
        for (j = 0; j < 3; j++) 
            for (k = 0; k < 3; k++) {
                FE[map[j]] += inverse[j][k] * tempE[j] / J;
                FU[map[j]] += inverse[j][k] * tempU[j] / J;
                FV[map[j]] += inverse[j][k] * tempV[j] / J;
            }
    }
}

void computeHeun(myStruct *tsunami, double *u, double *v, double *e, double *fu, double *fv, double *fe) {
    tsunami->U  = u;
    tsunami->V  = v;
    tsunami->E  = e;
    tsunami->FU = fu;
    tsunami->FV = fv;
    tsunami->FE = fe;
    
    integrateOnTriangles(tsunami);
    integrateOnEdges(tsunami);
    Phi_iPhiJOnEdgeAndInverse(tsunami);
}

void HeunMethod(myStruct *tsunami, HeunStruct *Heun) {
    double *USol  = tsunami->USol;
    double *VSol  = tsunami->VSol;
    double *ESol  = tsunami->ESol;

    int nElem = tsunami->nElem;
    int i;

    double *K1E     = Heun->K1E;
    double *K2E     = Heun->K2E;
    double *SumE    = Heun->SumE;
    double *K1U     = Heun->K1U;
    double *K2U     = Heun->K2U;
    double *SumU    = Heun->SumU;
    double *K1V     = Heun->K1V;
    double *K2V     = Heun->K2V;
    double *SumV    = Heun->SumV;
    double dt       = Heun->dt;

    for (i = 0; i < 3*nElem+1; i++) {
        K1E[i] = 0;
        K1U[i] = 0;
        K1V[i] = 0;
        K2E[i] = 0;
        K2U[i] = 0;
        K2V[i] = 0;
    }
    
    computeHeun(tsunami,USol,VSol,ESol,K1U,K1V,K1E);
    for (i = 0; i < nElem*3+1; i++) {
        SumE[i] = ESol[i] + dt * K1E[i];
        SumU[i] = USol[i] + dt * K1U[i];
        SumV[i] = VSol[i] + dt * K1V[i];
    }

    computeHeun(tsunami, SumU, SumV, SumE, K2U, K2V, K2E);
    for (i = 0; i < nElem*3+1; i++) {
        ESol[i] = ESol[i] + dt * (K1E[i] + K2E[i]) / 2.;
        USol[i] = USol[i] + dt * (K1U[i] + K2U[i]) / 2.;
        VSol[i] = VSol[i] + dt * (K1V[i] + K2V[i]) / 2.;
    }

}

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName) { 

    myStruct *tsunami = malloc(sizeof(myStruct));
    createStructureBasedOnFile(meshFileName, tsunami);
    
    int i,j;
    int nElem = tsunami->nElem;
    int *elem = tsunami->elem;
    
    double *E  = malloc(sizeof(double)*(3*nElem+1));
    for (i = 0; i < nElem; i++)
        for (j = 0; j < 3; j++)
            E[i*3+j] = tsunamiInitialConditionOkada(tsunami->X[elem[i*3+j]], tsunami->Y[elem[i*3+j]]);
    tsunami->ESol = E;

    HeunStruct *Heun = malloc(sizeof(HeunStruct));

    Heun->dt   = dt;
    Heun->K1E  = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->K2E  = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->SumE = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->K1U  = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->K2U  = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->SumU = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->K1V  = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->K2V  = malloc(sizeof(double) * (nElem*3 + 1));
    Heun->SumV = malloc(sizeof(double) * (nElem*3 + 1));

    for (i = 0; i <= nmax; i++) {
        HeunMethod(tsunami, Heun);
        if (i % sub == 0) tsunamiWriteFile(baseResultName, i, tsunami->USol, tsunami->VSol, tsunami->ESol, nElem, 3);
    }
        
    freeEverything(tsunami, Heun);
}
