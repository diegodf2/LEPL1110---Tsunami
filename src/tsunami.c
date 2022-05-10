# include "tsunami.h"

double tsunamiInitialConditionOkada(double x, double y)
{
  	double x3d = 4*R*R*x / (4*R*R + x*x + y*y);
  	double y3d = 4*R*R*y / (4*R*R + x*x + y*y);
  	double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  	double lat = asin(z3d/R)*180/PI;
  	double lon = atan2(y3d,x3d)*180/PI;
  	double lonMin = 142;
    double lonMax = 143.75;
    double latMin = 35.9;
    double latMax = 39.5;
  	double olon = (lonMin+lonMax)/2;
  	double olat = (latMin+latMax)/2;
  	double angle = -12.95*PI/180; 
  	double lon2 = olon + (lon-olon)*cos(angle) + (lat-olat)*sin(angle);
  	double lat2 = olat - (lon-olon)*sin(angle) + (lat-olat)*cos(angle);
  	if ( lon2 <= lonMax && lon2 >= lonMin && 
         lat2 >= latMin && lat2 <= latMax )	
    		return 1.0;
  	else	return 0.0;
}

void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub)
{
    int i,j;
    const char *basename = "%s-%08d.txt";
    
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    printf("%s\n",filename);
    FILE* file = fopen(filename,"w");
    
    fprintf(file, "Number of elem %d \n", nelem);
    fprintf(file, "Number of local values per element %d \n", nsub);
    for (i = 0; i < nelem; ++i) {
        for (j = 0; j < nsub; ++j) {
            int index = i*nsub+j;
            fprintf(file,"%d;%d;%le;%le;%le;\n",i,j,U[index],V[index],E[index]); }}    
    fclose(file);
}

void tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem)
{
    int i,j,trash,nelemFile,nsub;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"r");   
         
    fscanf(file, "Number of elem %d \n", &nelemFile);
    fscanf(file, "Number of local values per element %d \n", &nsub);     
    if (nelem != nelemFile) {
        printf("Error : wrong data file %d %d:-) \n",nelem,nelemFile);
        exit(0); }
    for (i = 0; i < nelem; ++i) {
    	  for (j = 0; j < nsub; ++j) {
        	  int index = i*nsub+j;
        	  fscanf(file,"%d;%d;%le;%le;%le;\n",&trash,&trash,&U[index],&V[index],&E[index]); }}                         
    fclose(file);
}

void tsunamiAnimate(double dt, int nmax, int sub,  const char *meshFileName, const char *baseResultName)
{
    int nElem,nNode,i,j,index,trash,*elem;
    double *X,*Y,*H,*E,*U,*V;
    int width,height;
    double mouse;
    
    int order = 1;
    int nout = sub;
    double t;
    double R = 6371220;
    double BathMax = 9368;
    GLfloat colors[9], coord[9];
   
    FILE* file = fopen(meshFileName,"r");
    if (file == NULL) {
    	  printf("Error : cannot open mesh file :-) \n");
        exit(0); }
 	  fscanf(file, "Number of nodes %d \n", &nNode);
  	X = malloc(sizeof(double)*nNode);
  	Y = malloc(sizeof(double)*nNode);
  	H = malloc(sizeof(double)*nNode);
	  for (i = 0; i < nNode; i++) 
    	  fscanf(file,"%d : %le %le %le  \n",&trash,&X[i],&Y[i],&H[i]); 
    fscanf(file, "Number of triangles %d \n", &nElem); 
  	elem = malloc(sizeof(int)*3*nElem);
 	  U = malloc(sizeof(double)*3*nElem);
 	  V = malloc(sizeof(double)*3*nElem);
 	  E = malloc(sizeof(double)*3*nElem);
  	for (i = 0; i < nElem; i++)     
    	  fscanf(file,"%d : %d %d %d \n", &trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]); 
  	fclose(file);
    

    glfwInit();
    GLFWwindow* window = glfwCreateWindow(1280,1280,"LEPL1110 : Tsunami project :-)",NULL,NULL);    
    glfwMakeContextCurrent(window);
    glShadeModel(GL_SMOOTH);
    glfwSwapInterval(1);

    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 8.0, 8.0, 8.0, 0.0 };
    
	  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    GLfloat light_radiance[] = {1., 1., 1., 1.};

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);	

    double t0 = 0;
    int frame0,frame = -1;
  
    do {
        
        t = glfwGetTime();                  
        frame0 = frame;
    	  frame = (int) ((t-t0) * 2);    
        if (frame0 != frame) {
        
            glfwGetCursorPos(window, &mouse, NULL );    mouse = 389;  
                
            const char *basename = "%s-%08d.txt";
            char filename[256];
            sprintf(filename,basename,baseResultName,frame*nout);
           	if (access(filename, F_OK) == 0) {
                printf("===  Reading local file %s %d %f \n",filename,frame,t);
                tsunamiReadFile(baseResultName,frame*nout,U,V,E,nElem); }
            
            glfwGetFramebufferSize(window,&width,&height);
            height = height > 0 ? height : 1;
            glViewport( 0, 0, width, height );

            glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);

            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            gluLookAt(0.0f,1.0f,0.0f,0.0f, 20.0f, 0.0f,0.0f,0.0f,1.0f);  
            glTranslatef(0.0f,14.0f,0.0f);
            double tt = 0;
            glRotatef(0.3f*(GLfloat)mouse + (GLfloat)tt*10.0f,0.0f,0.0f,1.0f);
            
            GLUquadricObj *quadratic = gluNewQuadric();         
            gluQuadricNormals(quadratic, GLU_SMOOTH); 
            glColor3f(1.0,1.0,1.0);
            gluSphere(quadratic,5.95,400,200);
            
            // A commenter pour supprimer la sphere interieure
            // Conseille pour les maillages grossiers :-)
     
            for (i=0; i < nElem; ++i) {
                for (j=0; j < 3; ++j) {
                    index = elem[3*i+j];
                    double value = H[index]/BathMax;
                    
            // Commenter la ligne qui suit pour visualiser la bathym�trie :-)
            // Puisque ligne qui pr�c�de contient la bathym�trie..
                             
                    value = E[3*i+j]*10;
                    if (value < 0) value = 0;
                    if (value > 1) value = 1; 
                    colors[j*3+0] = 3.5*(value)*(value);
                    colors[j*3+1] = (1-value)*(value)*3.5;
                    colors[j*3+2] = (1-value)*(1-value);
                    double x = X[index]; 
                    double y = Y[index]; 
                    double Factor = (4*R*R + x*x + y*y)*(R/6);
                    coord[j*3+0] = 4*R*R * x / Factor;
                    coord[j*3+1] = 4*R*R * y / Factor;
                    coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;  } 
      
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                glDisableClientState(GL_NORMAL_ARRAY);    
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);      
                                 
                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                for (j=0; j < 9; ++j)
                     coord[j] = coord[j] * 1.001;
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glDrawArrays(GL_LINE_LOOP, 0, 3);
                glDisableClientState(GL_VERTEX_ARRAY); }
        glfwSwapBuffers(window);
        glfwPollEvents(); }}
    while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
    
    glfwTerminate();
    exit( EXIT_SUCCESS );
}
