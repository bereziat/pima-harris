#include <Imlib2.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// pour compiler : gcc moravec.c -o moravec `imlib2-config --cflags --libs` -lm

////////recupération et création d'image

unsigned char * read_grayscale(char *fname, int *dimx, int *dimy)
{
  Imlib_Image *image;

  image = (Imlib_Image *)imlib_load_image(fname);
  if(image)
    {
      unsigned char *buf, *data;
      int l,c;

      imlib_context_set_image( image);
      *dimx = imlib_image_get_width();
      *dimy = imlib_image_get_height();
      data  = (unsigned char *)imlib_image_get_data();
      buf = (unsigned char*) malloc(sizeof(char)**dimx**dimy);
      for(l=0;l<*dimy;l++)
	for(c=0;c<*dimx;c++)
	  {
	    *buf++ = (*data + *(data+1) + *(data+2))/3;
	    data += 4;
	  }
      imlib_free_image();
      return buf-*dimx**dimy;
    }
  else
    {
      fprintf( stderr, "Erreur: image %s n'a pu etre ouverte.\n", fname);
      exit(33);
    }
}

double* read_grayscale_double(char *fname, int *dimx, int *dimy){
  unsigned char *imint = read_grayscale(fname, dimx, dimy);
  int tx=*dimx;
  int ty=*dimy;
  double* imdouble = (double *)malloc(sizeof(double)* tx* ty);
  int i;
  for(i=0;i<tx * ty; i++){
    imdouble[i]=((double)imint[i]) /255;
  }
  free(imint);
  return imdouble;
}


void write_grayscale(char *fname, int dimx, int dimy, unsigned char *buf)
{
  Imlib_Image *image;
  unsigned char *data;
  int l, c;
  char *ext;

  if( !(ext = strchr( fname, '.')))
    {
      fprintf( stderr, "Erreur: format image non reconnu\n");
      exit(34);
    }

  data = (unsigned char *)malloc(sizeof(char)*dimx*dimy*4);
  for( l=0; l<dimy; l++)
    for( c=0; c<dimx; c++)
      {
	*data = *buf;
	*(data+1)= *buf;
	*(data+2) = *buf;
	*(data+3) = 0;
	buf ++; data += 4;
      }
  data -= dimx*dimy*4;
  image = (Imlib_Image *)imlib_create_image_using_data( dimx, dimy, (DATA32*)data);
  imlib_context_set_image( image);
  imlib_image_set_format(ext+1);
  imlib_save_image( fname);
  imlib_free_image();
  free( data);
}
void write_grayscale_double(char *fname, int dimx, int dimy, double *imdouble){
  int * imint = malloc(sizeof(int)*dimx*dimy);
  unsigned char* imuchar = malloc(sizeof(unsigned char)*dimx*dimy);
  int i;
  for(i=0;i<dimx * dimy; i++){
    imint[i]=(int)(imdouble[i]*255);
    if(imint[i]<0){
      imint[i]=-imint[i];
    }
    if(imint[i]>255){
      imint[i]=255;
    }
    imuchar[i]=imint[i];
  }
  write_grayscale(fname, dimx, dimy, imuchar);
  free(imint);
  free(imuchar);
}

void write_color(char *fname, int dimx, int dimy, unsigned char *R, unsigned char *G, unsigned char *B) {
  Imlib_Image *image;
  unsigned char *data;
  char *ext;
  int l;
  
  if( !(ext = strchr( fname, '.'))) {
    fprintf( stderr, "Erreur: format image non reconnu\n");
    exit(34);
  }
  
  data = (unsigned char *)malloc(sizeof(char)*dimx*dimy*4);
  if( !data) {
    fprintf( stderr, "Erreur: plus de memoire libre (allocation de %ld octets).\n",
	     dimx*dimy*sizeof(float));
    exit(36);
  }
  for( l=0; l<dimy; l++) {
    int c;
    for( c=0; c<dimx; c++) {
      *data = *R++;
      *(data+1)= *G++;
      *(data+2) = *B++;
      *(data+3) = 0;
      data += 4;
    }
  }
  data -= dimx*dimy*4;
  image = imlib_create_image_using_data( dimx, dimy, (DATA32*)data);
  imlib_context_set_image( image);
  imlib_image_set_format( ext+1);
  imlib_save_image( fname);
  imlib_free_image();
  free( data);  
}

void write_grayscale_double_withe_color(char *fname, int dimx, int dimy, double *imdouble, double *imdoublecolore, double seuil){
  int * imint = malloc(sizeof(int)*dimx*dimy);

  unsigned char* imucharR = malloc(sizeof(unsigned char)*dimx*dimy);
  unsigned char* imucharG = malloc(sizeof(unsigned char)*dimx*dimy);
  unsigned char* imucharB = malloc(sizeof(unsigned char)*dimx*dimy);
  int i;
  for(i=0;i<dimx * dimy; i++){
		
    //imintR[i]=(int)((imdouble[i]+imdoublecolore[i])*255);
    imint[i]=(int)(imdouble[i]*255);
    if(imint[i]<0){
      imint[i]=-imint[i];
    }
    if(imint[i]>255){
      imint[i]=255;
    }

    imucharR[i]=imint[i];
    imucharG[i]=imint[i];
    imucharB[i]=imint[i];
		
    if(imdoublecolore[i]>seuil){
      imucharR[i]=255;
      imucharG[i]=0;
      imucharB[i]=0;
    }
  }

	
  write_color(fname, dimx, dimy, imucharB, imucharG, imucharR);
  free(imint);
  free(imucharR);
  free(imucharG);
  free(imucharB);
}

void write_grayscale_double_marque(char *fname, int dimx, int dimy, double *imdouble, int *imcolore){
	
	
  int * imint = malloc(sizeof(int)*dimx*dimy);

  unsigned char* imucharR = malloc(sizeof(unsigned char)*dimx*dimy);
  unsigned char* imucharG = malloc(sizeof(unsigned char)*dimx*dimy);
  unsigned char* imucharB = malloc(sizeof(unsigned char)*dimx*dimy);
  int i;
  for(i=0;i<dimx * dimy; i++){
		
    //imintR[i]=(int)((imdouble[i]+imdoublecolore[i])*255);
    imint[i]=(int)(imdouble[i]*255);
    if(imint[i]<0){
      imint[i]=-imint[i];
    }
    if(imint[i]>255){
      imint[i]=255;
    }

    imucharR[i]=imint[i];
    imucharG[i]=imint[i];
    imucharB[i]=imint[i];
		
    if(imcolore[i]==1){
      imucharR[i]=255;
      imucharG[i]=0;
      imucharB[i]=0;
    }
    if(imcolore[i]==2){
      imucharR[i]=0;
      imucharG[i]=255;
      imucharB[i]=0;
    }
    if(imcolore[i]==3){
      imucharR[i]=0;
      imucharG[i]=0;
      imucharB[i]=255;
    }
		
		
  }

	
  write_color(fname, dimx, dimy, imucharB, imucharG, imucharR);
  free(imint);
  free(imucharR);
  free(imucharG);
  free(imucharB);
}

///////// filtre et convolution

double* filtreGauss(unsigned int diametre, double sigma) {
  if(diametre%2 != 1){ 
    printf("erreur de taille de filtre\n");
    return NULL;
  }
  int i,u,v;
  double* tab_filtre = malloc(sizeof(double)*diametre*diametre);
  double somme=0;
  for(i=0; i<diametre*diametre; i++) {
    u = i%diametre - (diametre/2); 
    v = i/diametre - (diametre/2);
    tab_filtre[i] = exp((double)((u*u +v*v)/(2*sigma*sigma)));
    somme+=tab_filtre[i];
  }
  for(i=0; i<diametre*diametre; i++) {
    tab_filtre[i]=tab_filtre[i]/somme;
  }
  return tab_filtre;
}

double* filtreGT(unsigned int diametre, double amplification) {
  //gros rond
  if(diametre%2 != 1){ 
    printf("erreur de taille de filtre\n");
    return NULL;
  }
  double* tab_filtre = malloc(sizeof(double)*diametre*diametre);
  int i,u,v;
  for(i=0; i<diametre*diametre; i++) {
    u = i%diametre - (diametre/2);
    v = i/diametre - (diametre/2);
    if(u*u+v*v<=diametre*diametre){
      tab_filtre[i] = amplification;
    }else{
			
      tab_filtre[i] =0;
    }
  }
  return tab_filtre;
}



double* convolution(double* im, int dimx, int dimy, double* filtre, int diametrefiltre){
  int mf = diametrefiltre/2; //millieu filtre
  int i,j, k,l;
  double* res = calloc(dimx*dimy, sizeof(double));
  //ignor les bords, et les laisse à nul
  for(i=mf; i<dimx-mf; i++){
    for(j=mf; j<dimy-mf; j++){
      for(k=-mf;k<=mf;k++){
	for(l=-mf;l<=mf;l++){
	  res[i+j*dimx] += (im[(i+k)+(j+l)*dimx] * filtre[k+l*diametrefiltre]) ;
	}
      }
    }
  }
  return res ;
}


double* convolutionGauss(double* im, int dimx, int dimy){
  double* filtre = filtreGauss(7 , 2);
  double* res=convolution(im, dimx, dimy, filtre, 7);
  free(filtre);
  return res;
}


///////étapes de traitemant
 

void creetabXY(double* X, double* Y, double* image, int tl, int tc){
  //tl : taille des ligne ; tc : taille des colones
  int nx, ny;
  for(ny=0;ny<tc;ny++){
    X[ny*tl]= (image[ny*tl+1]-image[ny*tl])*2;
    for(nx=1;nx<tl-1;nx++){
      X[ny*tl+nx]=image[ny*tl+nx+1]-image[ny*tl+nx-1];
    }
    X[ny*tl+tl]= (image[ny*tl+tl]-image[ny*tl+tl-1])*2;
  }
  for(nx=1;nx<tl-1;nx++){
    Y[nx]=(image[nx+tl]-image[nx])*2;
    for(ny=1;ny<tc-1;ny++){
      Y[ny*tl+nx]=(image[ny*tl+nx-tl]-image[ny*tl+nx+tl]);
    }
    Y[tc*tl-tl+nx]=(image[tc*tl-tl+nx]-image[tc*tl-tl+nx-tl])*2;
  }
}

void creetabABC(double* A, double* B, double* C, double* X, double* Y, int dimx, int dimy){
  int i;
  for(i=0;i<dimx*dimy ; i++){
    A[i]=X[i]*X[i];
    B[i]=Y[i]*Y[i];
    C[i]=X[i]*Y[i];
  }
  int df = 3*3;
  double* filtre = filtreGauss(df , 3);
  double*tma = convolution(A, dimx, dimy, filtre, df);
  double*tmb = convolution(B, dimx, dimy, filtre, df);
  double*tmc = convolution(C, dimx, dimy, filtre, df);
	
  for(i=0;i<dimx*dimy ; i++){
    A[i]=tma[i];
    B[i]=tmb[i];
    C[i]=tmc[i];
  }
  free(tma);
  free(tmb);
  free(tmc);
}

void creetabAlphaBeta(double* Alpha, double* Beta, double* A, double* B, double* C, int dimx, int dimy){
  int i;
  float ab;
  float racine;
  for(i=0;i<dimx*dimy ; i++){
    ab=(A[i]+B[i])/2;
    racine=sqrt((A[i]-B[i])*(A[i]-B[i]) - C[i]*C[i])/2;
    Alpha[i]=ab-racine;
    Beta[i]=ab+racine;
  }
	

	
}

void creetabCoinBord(int* repereCB,double* Alpha, double* Beta,int dimx, int dimy, double seuilbord, double seuilplat){
  int i;
  double r;
  for(i=0;i<dimx*dimy ; i++){
    r= Alpha[i]*Beta[i] - seuilbord* (Alpha[i] + Beta[i] )*(Alpha[i] + Beta[i] );
    if(r<seuilplat && r>-seuilplat){
      repereCB[i]=0;
    }else{
      if(r>0){
	repereCB[i]=2;
      }else{
	repereCB[i]=3;
      }
    }
  }
  int k,l;
  for(i=0;i<dimx*dimy ; i++){
    if(repereCB[i]==2 && i>2*dimx && i<(dimy-2)*dimx && i%dimx >= 2 && i%dimx < dimx-2){
      for(k=-2;k<=2;k++){
	for(l=-2;l<=2;l++){
	  if(repereCB[i+k*dimx+l] != 2){ ////// ajout
	    repereCB[i+k*dimx+l]=1;
	  }
	}
      }
    }
  }
}

////// main


int main(int argc, char **argv){
  char *inName; //, *outName;
  //double seuil;
  if(argc > 2){
    
    double *dbuf_in, *dbuf_outX, *dbuf_outY, *dbuf_outXY;
    int dimx, dimy; // nombre de colonnes, de lignes de l'image 
    
    int i;
    
    //parametres de la commande 
    inName=argv[1];
    
		
    //lecture image 
    dbuf_in = read_grayscale_double(inName, &dimx, &dimy);
    
    // traitement sur l'image 
    dbuf_outY = (double *)malloc(sizeof(double)* dimx*dimy);
    dbuf_outX = (double *)malloc(sizeof(double)*dimx*dimy);
    dbuf_outXY = (double *)malloc(sizeof(double)*dimx*dimy);
    
    
    creetabXY(dbuf_outX,dbuf_outY,dbuf_in,dimx, dimy);
    for(i=0;i<dimx * dimy; i++){
      dbuf_outXY[i]=fabs(dbuf_outX[i])+fabs(dbuf_outY[i]);
    }
    
    double *A = (double *)malloc(sizeof(double)* dimx*dimy);
    double *B = (double *)malloc(sizeof(double)*dimx*dimy);
    double *C = (double *)malloc(sizeof(double)*dimx*dimy);
    
    creetabABC(A, B, C, dbuf_outX, dbuf_outY, dimx, dimy);
		
    double *Alpha = (double *)malloc(sizeof(double)*dimx*dimy);
    double *Beta = (double *)malloc(sizeof(double)*dimx*dimy);
		
    creetabAlphaBeta(Alpha, Beta, A, B, C, dimx, dimy);
		
		
    double * fGT = filtreGT(9,15);
    double *Alphaample = convolution(Alpha, dimx, dimy, fGT, 9);
		
    int* repereCB = malloc(sizeof(int)*dimx*dimy);
    creetabCoinBord( repereCB,Alpha,Beta, dimx,dimy, 0.04 , atof(argv[2]));
		
    // ecriture image 
		
    /*  // image des tableau intermédiaire
	write_grayscale_double("tabX.jpg", dimx, dimy, dbuf_outX);
	write_grayscale_double("tabY.jpg", dimx, dimy, dbuf_outY);
	write_grayscale_double("tabXY.jpg", dimx, dimy, dbuf_outXY);
		
	write_grayscale_double("tabA.jpg", dimx, dimy, A);
	write_grayscale_double("tabB.jpg", dimx, dimy, B);
	write_grayscale_double("tabC.jpg", dimx, dimy, C);
		
	write_grayscale_double("tabAlpha.jpg", dimx, dimy, Alpha);
	write_grayscale_double("tabBeta.jpg", dimx, dimy, Beta);
	write_grayscale_double("tabAlphaample.jpg", dimx, dimy, Alphaample);
		
	write_grayscale_double_withe_color("taborigine.jpg", dimx, dimy, dbuf_in, dbuf_outXY, 0.2 );
    */
		
    write_grayscale_double_marque("imageFinal.jpg", dimx, dimy, dbuf_in,repereCB);
    free(dbuf_in);
    free(dbuf_outX);
    free(dbuf_outY);
    free(dbuf_outXY);
    free(A);
    free(B);
    free(C);
		
    free(Alpha);
    free(Beta);
		
    free(fGT);
    free(Alphaample);
    free(repereCB);

  }else{
    printf("Usage: %s image-in\n", *argv);
  }
	
  return 0;
}
