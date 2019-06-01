//NAME:: Vis_Avg.c: Does Fourier transform ( DFT) of an image and creates visibilities at desired baselines from a uvfits file
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include <sys/stat.h>
# include <nr.h>
# include <nrutil.h>

# define DR (M_PI/180.)

// INCODE parameters
fftw_complex *GVis;

int    NGridI;
float  DIx, DIy;
float BMaj, BMin;

int Access(char *filename){
  struct stat   buffer;   
  return (stat (filename, &buffer) == 0);
}

void printerror(int status){
  if (status){
    fits_report_error(stderr, status); // Print error report
    exit( status );    // Terminate the program, returning error status
  }
}
int Read_Image_FITS(char *INFILE){
  
  fitsfile *fptr;  
  long ii, jj, index, index1;
  int  status,  anynull, nfound;
  long fpixel, npixels, *naxes, naxis;
  double  *dataf, *Idata, tx, ty, nullval; 
  fftw_plan dft_r2i;
  char comment[FLEN_KEYWORD];
 
  printf("Reading Image FITS file...\n");
  status = 0;
  if(fits_open_file(&fptr, INFILE, READONLY, &status))
    return status;
  if(fits_read_key_lng(fptr, "NAXIS", &naxis, comment, &status))
    return status;
  naxes = (long *)calloc(naxis, sizeof(long));
  if(fits_read_keys_lng(fptr, "NAXIS", 1, naxis, naxes, &nfound, &status))
    return status;
  
 if(naxes[0]!=naxes[1]){
    printf("\nWARNING: Image is not symmetric in x and y [%ld\t%ld]\n\n", naxes[0], naxes[1]); 
    return EXIT_FAILURE;
  }
  NGridI = naxes[0];
  printf("\nNGridI = %d\n", NGridI);

  
  if(fits_read_key_flt(fptr, "CDELT1", &DIx, comment, &status)){
    printf("Error reading CDELT1\n");
    return status;
  }
  if(fits_read_key_flt(fptr, "CDELT2", &DIy, comment, &status)){
    printf("Error reading CDELT2\n");
    return status;
  }
 if(fits_read_key_flt(fptr, "BMAJ", &BMaj, comment, &status)){
    printf("Error reading BMaj\n");
    return status;
  }
 if(fits_read_key_flt(fptr, "BMIN", &BMin, comment, &status)){
    printf("Error reading BMin\n");
    return status;
  }
 
  printf("DIx    = %.2e arc sec\t", DIx*3600.);
  printf("DIy    = %.2e arc sec\t", DIy*3600.);

  printf("BMaj   = %.2e arc sec\t", BMaj*3600.);
  printf("BMin   = %.2e arc sec\t", BMin*3600.);
  
  DIx  *= DR;
  DIy  *= DR;
  BMaj *= DR;
  BMin *= DR;
  
  //// here comes new things
  nullval  = 0;
  npixels  = NGridI*NGridI;
  fpixel = 1;
  
  Idata = (double *)calloc(npixels, sizeof(double));
  dataf = (double *)calloc(NGridI*(NGridI+2), sizeof(double));
  GVis = (fftw_complex *)&dataf[0];
  dft_r2i = fftw_plan_dft_r2c_2d (NGridI, NGridI, dataf, GVis, FFTW_ESTIMATE);

  if(fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, Idata, &anynull, &status))
    return status;
  double Norm;
  
  Norm = 1./sqrt(1.*NGridI*NGridI);

  // Normalization below consider the followings:
  // BMAJ = 2 x sqrt ( 2 ln (2) ) sigma_MAJ etc
  // Area of the beam = 2 * pi * sigma_MAJ * sigma_MIN
  // Hence: Area = 2 * pi * BMAJ *BMIN  / (2 / sqrt  ( 2 * log (2)) )^2
  //             = pi * BMAJ * BMIN / ( 4 * log(2))
  
  Norm = Norm*(M_PI*BMaj*BMin/DIx/DIy/4./log(2.));

  printf("Norm =  %e\n", Norm);

  for(jj=0;jj<NGridI;++jj)
    for(ii=0;ii<NGridI;++ii){
      
      index=jj*NGridI+ii;
      index1=ii*(NGridI+2)+jj;
      dataf[index1] = Idata[index];
    }
        
  fftw_execute(dft_r2i);
    
  for(ii=0;ii<NGridI;++ii)
    for(jj=0;jj<=NGridI/2;++jj){
	
      index1 = ii*(NGridI/2+1)+jj;
      index  = ii*(NGridI/2+1)+jj;
      // Following lines does fftshift
      GVis[index1][0] = Norm*pow(-1.,ii+jj)*GVis[index][0];
      GVis[index1][1] = Norm*pow(-1.,ii+jj)*(-1.*GVis[index][1]);
    }
  
  fftw_destroy_plan(dft_r2i);
  if(fits_close_file(fptr, &status))
    return status;
  
  free(Idata);
  
  return EXIT_SUCCESS;
}
int main(int argc, char *argv[]){

  int STATFLAG;
  char INVISFITS[128], INIMFITS[128], OUTFITS[128];
 
  if(argc!=4){

      printf("Usage: %s <input UVFITS file> <input IMAGE fits file> <Output Fits File> \n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INVISFITS);
  sscanf(argv[2],"%s",INIMFITS);
  sscanf(argv[3],"%s",OUTFITS);

  BMaj = DR*BMaj;
  BMin = DR*BMin;

  int BChan, EChan, ii_G;
  long  OChan;

  if(Access(INVISFITS)==0){
    printf("Input File %s does not exists\n", INVISFITS);
    return EXIT_FAILURE;
  }
  if(Access(INIMFITS)==0){
    printf("Input File %s does not exists\n", INIMFITS);
    return EXIT_FAILURE;
  }
  if(Access(OUTFITS)!=0){
    printf("Output File %s exists\n", OUTFITS);
    return EXIT_FAILURE;
  }
  Read_Image_FITS(INIMFITS);

    
  fitsfile *fptr,*fptr1, *fptro;
  int status=0,anynul=0;
  int ii, jj;
  float nulval=0.;
  long GCount, PCount, NStokes, NChan, NCmplx, group, nel;
  float del_chan, ii_r, nu_r, nu_G, *randpar;
  float *vel, *vis, *vmom0, *vmom1;
  char comment[FLEN_COMMENT];
  int stokes,chan,index, indexm;
  float Umax, Uval, dU, data[3];
  float uu, vv, uus, vvs, signv;
     
  fits_open_file(&fptr, INVISFITS, READONLY, &status);
  fits_create_file(&fptro, OUTFITS, &status);
  fits_copy_header(fptr, fptro, &status);

  if(fits_read_key_lng(fptr, "NAXIS2", &NCmplx, comment, &status))
    printerror( status );
  printf("NAXIS2=%ld\n", NCmplx);
  
  if(fits_read_key_lng(fptr,"NAXIS4", &NStokes,comment,&status))
    printerror( status );
  printf("NAXIS3=%ld\n", NStokes);
  
  if(fits_read_key_lng(fptr,"NAXIS5", &NChan, comment,&status))
    printerror( status );
  printf("NAXIS4=%ld\n", NChan);
  
  if(fits_read_key_lng(fptr,"GCOUNT", &GCount,comment,&status))
    printerror( status );
  printf("GCOUNT=%ld\n", GCount);
  
  if(fits_read_key_lng(fptr,"PCOUNT", &PCount, comment,&status))
    printerror( status );
  printf("PCOUNT=%ld\n", PCount);
 
  nel=NCmplx*NStokes;
  printf("nel=%ld\n",nel);

  if(NStokes!=1){
    fprintf(stderr, "Dont work with multi-stokes files\n");
    return EXIT_FAILURE;
  }
  if(NChan!=1){
    fprintf(stderr, "Dont work with multi-channel files\n");
    return EXIT_FAILURE;
  }

  randpar = (float*)calloc(PCount, sizeof(float));
  float DL; 
  ////////////read header complete/////////////////
  DL = (fabs(DIx) <= fabs(DIy)) ? fabs(DIx) : fabs(DIy);
  dU   = 1./(1.*NGridI*DL);
  Umax = 1./(2.*DL);
  printf("Umax = %e\n", Umax);
  for(group=1;group<=GCount;group++){
      
    if(fits_read_grppar_flt(fptr, group, 1, PCount, randpar, &status ) )
      printerror(status);
    if(fits_write_grppar_flt(fptro, group, 1, PCount, randpar, &status ) )
      printerror(status);

    uu = randpar[0];
    vv = randpar[1];

    Uval=sqrt(uu*uu+vv*vv);
    if(Uval<Umax){
    
    signv = (vv<0.) ? -1. : 1.;
    uus   = signv*uu;
    vvs   = signv*vv;
   
    ii  = (int)roundf(uus/dU);
    ii  = (ii <0) ? NGridI+ii : ii;
    jj  = (int)roundf(vvs/dU);
    
    index = (NGridI/2+1)*ii+jj;
    
    data[0] = GVis[index][0];
    data[1] = (signv*GVis[index][1]);
    data[2] = 1.;
    }
    else{
      data[1] = 0.;
      data[2] = 0.;
      data[0] = 0.;
    }

    if(fits_write_img_flt(fptro, group, 1, 3, data, &status))
      printerror(status);
  }
  printf("Output written in %s\n", OUTFITS);
  fits_close_file(fptr,&status);
  fits_close_file(fptro,&status);

  return EXIT_SUCCESS;
}

