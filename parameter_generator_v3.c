#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>

void SimpleRandomNumberInitilaization(gsl_rng *,long int);
/*void OLDDomCycleReactIndex(int , int *, int *,int , int);*/
void DomCycleReactIndex(int , int *, int *,int , int , int);
void WriteKineticsRates(char *, double *, double , double , double , double , double , double , double , double , int , int , int , int );
void ReadKineticsRates(char *, double *, double *, double *, double *, double *, double *, double *, int , int , double *, double *,int *, int *);

gsl_rng *r;
FILE *input, *output, *fp;

int main(int argc, char *argv[])
{
  int nA=atoi(argv[1])/*3*/;
  int nB=atoi(argv[2])/*3*/;
  int kSetGeneration=atoi(argv[3]); /*0: no generating k_set or k_SET; 1: generating k_set; 2: generating k_SET*/
  int DominantCycle=atoi(argv[4]); /*<0: no catalitic aid, 0: system died out; 1: A is dominant cyckle; 2: B is dominant cyckle; 3: A and B cyckles coexist; 4: something wrong*/
  long int InitRandNum=atoi(argv[5]);
  int reverse=atoi(argv[6]); /*0: no reverse, 1: reverse*/
  double hanyados=atof(argv[7]); /* revers rate= toward rate/hanyados*/
  int HetCat=atoi(argv[8]); /*0: no heterocatalysis; 1: HC: sub1: A1 or B1, sub2: X; 2: rHC: sub1: X, sub2: A1 or B1 */
  int R=atoi(argv[9]); /** the number of reaction where the revers reaction is occure **/
  int TModel=atoi(argv[10]);
  int JoinedCycles=atoi(argv[11]); /** 0: no joined cycles, 1: cycles joined at JpinesdPointA=JoinedPointB**/
  int JoinedPointA=atoi(argv[12]); /** 0: no joined cycles, >0: cycle A is joined to cycle B at this point **/
  int JoinedPointB=atoi(argv[13]); /** 0: no joined cycles, >0: cycle B is joined to cycle A at this point **/

	
  double *Param=NULL;
  int LKset, NReact, NumOfEq=0, Enzyme=0, InteractionT2=0, InteractionT1=0, reacT=-9, reacR=-9, PozRev, pointE, pointP;
  double **reactions=NULL,kon1,koff1,kon2,koff2,kcat,kcatr, V, v, p;
  
  register int i;
  
  if(TModel==0)
  {
    V=v=0.00;
  }
  else
  {
    V=1.00;
    v=0.02;
  }

  PozRev=R-1;
/** **//*printf("nA = %d\nnB = %d\nkSetGeneration = %d\nDominantCycle = %d\nInitRandNum = %ld\nreverse = %d\nhanyados = %f\nHetCat = %d\nR(%d) = %d\nTModel = %d\n",nA, nB, kSetGeneration, DominantCycle, InitRandNum, reverse, hanyados, HetCat, R ,PozRev, TModel);*/
/** **//*exit(101);*/

  NumOfEq=NumOfEq;
  Enzyme=Enzyme;
	kon1=koff1=kon2=koff2=kcat=kcatr=0.00;
  NReact=nA+nB;
  reactions=(double**)calloc(NReact, sizeof(double*));
  for(i=0;i<NReact;i++)
    reactions[i]=(double*)calloc(6, sizeof(double));
  LKset=(2*(nA+nB));
  Param=(double*)calloc((unsigned int)LKset,sizeof(double));
  
  if(kSetGeneration==1)
  {
		r=gsl_rng_alloc(gsl_rng_mt19937);
		SimpleRandomNumberInitilaization(r,InitRandNum);
/** **//*printf("A random number [0,1]: %f\n",gsl_rng_uniform(r));*/
/** **/  /*printf("sor= %d, nA= %d, nB= %d, InteractionT1=%d, InteractionT2=%d\n", sor, nA, nB, InteractionT1, InteractionT2);*/
		/*LKset=(2*(nA+nB));
		Param=(double*)calloc(LKset,sizeof(double));*/

/** **/ /* NO REVERSE REACTION, NO CATALYTIC AID */
		for(i=0;i<(nA+nB);i++)
		{
			Param[i]=round(100*((10.0*gsl_rng_uniform(r))+0.1))/100.0;
			if(Param[i]>10.0)
				Param[i]=10.0;
			else;
			if(Param[i]<0.1)
				Param[i]=0.1;
			else;
		
			Param[i+nA+nB]=0.00;			
		}

		WriteKineticsRates("k_set", &Param[0], kon1, kon2, kcat, koff1, koff2, kcatr, v,V,NumOfEq, Enzyme, nA, nB);   
/** **//*exit(215);*/
  }
  else if(kSetGeneration==2)
  {

	r=gsl_rng_alloc(gsl_rng_mt19937);
	SimpleRandomNumberInitilaization(r,InitRandNum);
/** **//*printf("A random number [0,1]: %f\n",gsl_rng_uniform(r));*/
/** **/  /*printf("sor= %d, nA= %d, nB= %d, InteractionT1=%d, InteractionT2=%d\n", sor, nA, nB, InteractionT1, InteractionT2);*/
	/*LKset=(2*(nA+nB));
	Param=(double*)calloc(LKset,sizeof(double));*/
		
	ReadKineticsRates("k_set", &Param[0], &kon1, &kon2, &kcat, &koff1, &koff2, &kcatr, nA, nB, &V, &v, &NumOfEq, &Enzyme);
/** **//*printf("reverse: %d\n", reverse);*/
    
    if(reverse!=0)
    {
      /*OLDDomCycleReactIndex(DominantCycle,&reacT, &reacR,nA, nB); */
			DomCycleReactIndex(DominantCycle,&reacT, &reacR,nA, nB,PozRev); 
      Param[reacR]=Param[reacT]/hanyados;
    }
    else {;}
/** **//*exit(256);*/
    
    if(nB>0)
    {
      if(DominantCycle==1)
        p=1.0/(double)nB; /* the probability of Interaction Type 2 on cycle A */
      else if(DominantCycle==2)
        p=1.0/(double)nA; /* the probability of Interaction Type 2 on Cycle B */
      else p=0.5; /* There is no dominant cycle, but the bash file creates a new k-set before running this file! */
	
      if(gsl_rng_uniform(r)>=p)
      {
        InteractionT1=1;
        InteractionT2=0;
        if(DominantCycle==1)
        {
          Enzyme=gsl_rng_uniform_int(r,nA)+1;
          NumOfEq=nA+gsl_rng_uniform_int(r,nB-1)+2;
        }
        else if(DominantCycle==2)
        {
          Enzyme=nA+gsl_rng_uniform_int(r,nB)+1;
          NumOfEq=gsl_rng_uniform_int(r,nA-1)+2;
        }
        else;
			
        kon1=round(100*((10.0*gsl_rng_uniform(r))+0.1))/100.0;
        koff1=kon1/((double)(2.0+(gsl_rng_uniform(r)*8.0)));
        kcat=round(100*((10.0*gsl_rng_uniform(r))+0.1))/100.0;
        kcatr=0.00;	
			
        WriteKineticsRates("k_SET", &Param[0], kon1,kon2,kcat, koff1, koff2, kcatr, v,V,NumOfEq, Enzyme, nA, nB); 	
			
      }
      else
      {
        InteractionT1=0;
        InteractionT2=1;
        if(DominantCycle==1)
        {
          Enzyme=gsl_rng_uniform_int(r,nA)+1;
          NumOfEq=nA+1;
        }
        else if(DominantCycle==2)
        {
          Enzyme=nA+gsl_rng_uniform_int(r,nB)+1;
          NumOfEq=1;
        }
        else;
			
        kon1=round(100*((10.0*gsl_rng_uniform(r))+0.1))/100.0;
        koff1=kon1/((double)(2.0+(gsl_rng_uniform(r)*8.0)));
        kon2=round(100*((10.0*gsl_rng_uniform(r))+0.1))/100.0;
        koff2=kon2/((double)(2.0+(gsl_rng_uniform(r)*8.0)));
        kcat=round(100*((10.0*gsl_rng_uniform(r))+0.1))/100.0;
        kcatr=0.00;
			
        WriteKineticsRates("k_SET", &Param[0], kon1,kon2,kcat, koff1, koff2, kcatr, v,V,NumOfEq, Enzyme, nA, nB); 
      }
    }
/** **//*exit(156);*/
  }
  else if(kSetGeneration==0)
  {
    /*LKset=(2*(nA+nB));*/
    /*Param=(double*)calloc(LKset,sizeof(double));*/

    if(HetCat==0)
    {
      ReadKineticsRates("k_set", &Param[0], &kon1, &kon2, &kcat, &koff1, &koff2, &kcatr, nA, nB, &V, &v, &NumOfEq, &Enzyme);

      if(reverse!=0)
      {
        /*OLDDomCycleReactIndex(DominantCycle,&reacT, &reacR,nA, nB);*/
				DomCycleReactIndex(DominantCycle,&reacT, &reacR,nA, nB,PozRev); 
        Param[reacR]=Param[reacT]/hanyados;
      }
      else;   
      WriteKineticsRates("k_SET2", &Param[0], kon1, kon2, kcat, koff1, koff2, kcatr, v,V,NumOfEq, Enzyme, nA, nB);
/** **//*exit(5552);*/
    }
    else if(HetCat>0)
    {
      ReadKineticsRates("k_SET", &Param[0], &kon1, &kon2, &kcat, &koff1, &koff2, &kcatr, nA, nB, &V, &v, &NumOfEq, &Enzyme);
      
       if(reverse!=0)
      {
        /*OLDDomCycleReactIndex(DominantCycle,&reacT, &reacR,nA, nB); */
				DomCycleReactIndex(DominantCycle,&reacT, &reacR,nA, nB, PozRev); 
        Param[reacR]=Param[reacT]/hanyados;
      }
      else;
      
      WriteKineticsRates("k_SET2", &Param[0], kon1, kon2, kcat, koff1, koff2, kcatr, v,V,NumOfEq, Enzyme, nA, nB);
    }
    else
    {
      printf("Error: No information about the reading of k_set or k_SET!");
      exit(101);
    }
  }
  else
  {
    printf("Error: No information about k_set generartion!");
    exit(102);
  }

  if(NumOfEq==0)
  {
	InteractionT1=0;
	InteractionT2=0;
  }
  else if((NumOfEq==1)||(NumOfEq==(nA+1)))
  {
	InteractionT1=0;
	InteractionT2=1;
  }
  else
  {
	InteractionT1=1;
	InteractionT2=0;
  }


	
  output=fopen("proba_parameters.txt","w");
  fprintf(output,"# general info: Start, End, DeltaT, Num. of reactions,Num. of Chem. Entities, Sum of cc,  V, v, Tmodel, JoinedCycles, Num. of Interactions Type 1, Num. of Interactions Type 2\n");
  fprintf(output,"0.00,100000.00,1.0,%d,%d,0.00,%f,%f,%d,%d,%d,%d\n", nA+nB, nA+nB+1, V, v, TModel, JoinedCycles, InteractionT1, InteractionT2);
/** **//*printf("0.00,1000000.00,1.0,%d,%d,0.00,%f,%f,%d,0,%d,%d\n", nA+nB, nA+nB+1, V, v, TModel, InteractionT1, InteractionT2);*/
  fflush(output);
/** **//*exit(158);*/

  fprintf(output,"# inflow: index of chemical entity, regime (IF Tmodel==1: 0: unlimited, -Value: limited, +Value: constant inflow)\n");
  if(TModel==0)
    fprintf(output,"%d,0.00\n",nA+nB+1);
  else if(TModel==1)
    fprintf(output,"%d,1.00\n",nA+nB+1);
  else
  {
    printf("Error:  Simulation regime (unlimited or chemostat) is unknown! (1)\n");
    exit(12);
  }
  fflush(output);
/** **//*exit(158);*/

  fprintf(output,"# initial cc: index of chemical entity, cc\n");
  if(TModel==0)
  {
    fprintf(output,"%d,1.00\n",1);
    if(nB==0)
    {
      for(i=1;i<(nA+nB);i++)
        fprintf(output,"%d,0.00\n",i+1);
    }
    else if(nB>0)
    {
      for(i=1;i<(nA);i++)
        fprintf(output,"%d,0.00\n",i+1);
      if(JoinedPointB==nA+1)
        fprintf(output,"%d,0.00\n",nA+1);
      else fprintf(output,"%d,1.00\n",nA+1);
      for(i=(nA+1);i<(nA+nB);i++)
        fprintf(output,"%d,0.00\n",i+1);
    }
    else{;}
    fprintf(output,"%d,1.00\n",nA+nB+1);
  }
  else if(TModel==1)
  {
    for(i=0;i<(nA+nB);i++)
    {
      if(JoinedPointB==(i+1))
        fprintf(output,"%d,0.00\n",i+1);
      else fprintf(output,"%d,0.01\n",i+1);
    }
  }
  else
  {
    printf("Error:  Simulation regime (unlimited or chemostat) is unknown! (2)\n");
    exit(12);
  }
  fflush(output);
/** **//*exit(158);*/
  
  fprintf(output,"# reactions: Educt1, Educt2, k_on, k_rev, Product1, Product2 (k_on<0 and/or k_off<0 -> random rata generation)\n");
  fprintf(output,"1,%d,%.2f,%.2f,2,0\n", nA+nB+1, Param[0], Param[nA+nB]);
  reactions[0][0]=1;
  reactions[0][1]=(double)nA+nB+1;
  reactions[0][2]=Param[0];
  reactions[0][3]=Param[nA+nB];
  reactions[0][4]=2.00;
  reactions[0][5]=0.00;

  for(i=2;i<nA;i++)
  {
      fprintf(output,"%d,0,%.2f,%.2f,%d,0\n", i, Param[i-1], Param[nA+nB+i-1], i+1);
      reactions[i-1][0]=(double)i;
      reactions[i-1][1]=0.0;
      reactions[i-1][2]=Param[i-1];
      reactions[i-1][3]=Param[nA+nB+i-1];
      reactions[i-1][4]=(double)i+1;
      reactions[i-1][5]=0.00;
  }
  fprintf(output,"%d,0,%.2f,%.2f,1,1\n", nA, Param[nA-1], Param[nA+nB+nA-1]);
  reactions[nA-1][0]=(double)nA;
  reactions[nA-1][1]=0.0;
  reactions[nA-1][2]=Param[nA-1];
  reactions[nA-1][3]=Param[nA+nB+nA-1];
  reactions[nA-1][4]=(double)1;
  reactions[nA-1][5]=(double)1;

  if(nB>0)
  {
      if((nA+1)==JoinedPointB)
        pointE=JoinedPointA;
      else pointE = nA+1;
      if((nA+2)==JoinedPointB)
        pointP=JoinedPointA;
      else pointP = nA+2;

      fprintf(output,"%d,%d,%.2f,%.2f,%d,0\n", pointE, nA+nB+1, Param[nA], Param[nA+nB+nA], pointP);
      reactions[nA][0]=(double)pointE;
      reactions[nA][1]=(double)nA+nB+1;                   
      reactions[nA][2]=Param[nA];
      reactions[nA][3]=Param[nA+nB+nA];
      reactions[nA][4]=(double)pointP;
      reactions[nA][5]=0.00;
  
 
      for(i=2;i<nB;i++)
      {
        if((nA+i)==JoinedPointB)
          pointE=JoinedPointA;
        else pointE = nA+i;
        if((nA+i+1)==JoinedPointB)
          pointP=JoinedPointA;
        else pointP = nA+i+1;
        fprintf(output,"%d,0,%.2f,%.2f,%d,0\n", pointE, Param[nA+i-1], Param[nA+nB+nA+i-1], pointP);
        reactions[nA+i-1][0]=(double)pointE;
        reactions[nA+i-1][1]=0.0;
        reactions[nA+i-1][2]=Param[nA+i-1];
        reactions[nA+i-1][3]=Param[nA+nB+nA+i-1];
        reactions[nA+i-1][4]=(double)pointP;
        reactions[nA+i-1][5]=0.00;
      }

      if((nA+nB)==JoinedPointB)
        pointE=JoinedPointA;
      else pointE = nA+nB;
      if((nA+1)==JoinedPointB)
        pointP=JoinedPointA;
      else pointP = nA+1;
      fprintf(output,"%d,0,%.2f,%.2f,%d,%d\n", pointE, Param[(nA+nB)-1], Param[nA+nB+nA+nB-1], pointP,pointP);
      reactions[(nA+nB)-1][0]=(double)pointE;
      reactions[(nA+nB)-1][1]=0.0;
      reactions[(nA+nB)-1][2]=Param[(nA+nB)-1];
      reactions[(nA+nB)-1][3]=Param[nA+nB+nA+nB-1];
      reactions[(nA+nB)-1][4]=(double)pointP;
      reactions[(nA+nB)-1][5]=(double)pointP;
  }

  
/** **/  /*for(i=0;i<NReact;i++)*/
/** **/  /*{*/
/** **/    /*  for(j = 0; j<(nA+nB-1);j++) */
/** **/      /* printf("%f,", reactions[i][j]);*/
/** **/     /*printf("%f\n", reactions[i][nA+nB-1]);*/
/** **/  /*}*/
/** **/ /*exit(16);*/
  fflush(output);
/** **//*exit(158);*/
  fprintf(output,"# Interaction type 1: E+S1=ES1->E+P1+P2: S1, E, kon,koff, kcat, kcatr, P1, P2\n");
/** **/ /*exit(17);*/
	
  if((InteractionT1>0) /*&& (ExitRunningThisKSet==0)*/ )
	{
    fprintf(output,"%d,%d,%.2f,%.2f,%.2f,%.2f,%d,%d\n", (int)reactions[NumOfEq-1][0], Enzyme, kon1, koff1, kcat, kcatr, (int)reactions[NumOfEq-1][4], (int)reactions[NumOfEq-1][5]);
/** **//*printf("%d,%d,%.2f,%.2f,%.2f,%.2f,%d,%d\n", (int)reactions[NumOfEq-1][0], Enzyme, kon1, koff1, kcat, kcatr, (int)reactions[NumOfEq-1][4], (int)reactions[NumOfEq-1][5]);*/
	}
  else; 
  fflush(output);
/** **//*exit(158);*/

  fprintf(output,"# Interaction type 2: E+S1=ES1, ES1+S2=ES1S2->E+P1+P2: S1, S2, E, kon1, koff1, kon2, koff2, kcat, kcatr, P1, P2\n");
/** **/ /*exit(17);*/
  if((InteractionT2>0) /*&& (ExitRunningThisKSet==0)*/)
	{
    if(HetCat==1)
    {
      fprintf(output,"%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d\n",(int)reactions[NumOfEq-1][0], (int)reactions[NumOfEq-1][1], Enzyme, kon1, koff1, kon2, koff2, kcat, kcatr, (int)reactions[NumOfEq-1][4], (int)reactions[NumOfEq-1][5]);
    }
    else if(HetCat==2)
    {
      fprintf(output,"%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d\n",(int)reactions[NumOfEq-1][1], (int)reactions[NumOfEq-1][0], Enzyme, kon1, koff1, kon2, koff2, kcat, kcatr, (int)reactions[NumOfEq-1][4], (int)reactions[NumOfEq-1][5]);
    }
    else
    {
      printf("Error: No information about the order of subsrtrates order in heterocatalysis!");	
      exit(103);
    }
	}
  else;
  fflush(output);
/** **//*exit(158);*/
  fprintf(output,"# Joined cycles at these points\n");
  fprintf(output,"%d=%d\n", JoinedPointA, JoinedPointB);
  fflush(output);
  fclose(output);
	
	
  
  free(Param);
  for(i=0;i<NReact;i++)
    free(reactions[i]);
  free(reactions);
  gsl_rng_free (r);
  
  return 0;
}


void SimpleRandomNumberInitilaization(gsl_rng *r,long int InitVOfRN)
{
	long int dummy=0,RandNum=0;
/** **/ /*printf("### START SimpleRandomNumberInitilaization()\n");*/
/** **/	/*r=gsl_rng_alloc(gsl_rng_mt19937);*/
	
/** **/ /*printf("InitVOfRN= %ld\n",InitVOfRN);*/
	if(InitVOfRN==0)
	{
		fp = fopen("/dev/urandom", "r"); 
		dummy=fread(&RandNum, sizeof(long),1, fp);
		fclose(fp);
		dummy=dummy;
		if(RandNum<0)
			RandNum*=(-1);
		else;
/** **/ /*printf("RandNum= %ld\n",RandNum);*/
		gsl_rng_set(r,RandNum);
	}
	else
	{
		if(InitVOfRN<0)
			InitVOfRN*=(-1);
		else;
/** **/ /*printf("InitVOfRN= %ld\n",InitVOfRN);*/
		gsl_rng_set(r,InitVOfRN);
	}
	/** **/ /*printf("InitVOfRN= %ld, InitialRandNumber: %ld\n",InitVOfRN,RandNum);*/
	/** **/ /*printf("### END SimpleRandomNumberInitilaization()\n");*/
}

void   WriteKineticsRates(char *fileName, double *P, double kkon1, double kkon2, double kkcat, double kkoff1, double kkoff2, double kkcatr, double vv, double VV, int NOEq, int E, int nnA, int nnB)
{
    int i;
    
    output=fopen(fileName,"w");
    for(i=0;i<nnA+nnB;i++) /** REACTION RATES TOWARD**/
    {
      fprintf(output,"%.2f ",P[i]);
/** **/		/*printf("%.2f ",P[i]);*/
    }
    fprintf(output,"%.2f %.2f %.2f ", kkon1, kkon2, kkcat);
/** **/	/*printf("%.2f %.2f %.2f ", kkon1, kkon2, kkcat);*/
    for(i=0;i<nnA+nnB;i++) /** REACTION RATES BACKWARD**/
    {
      fprintf(output,"%.2f ",P[nnA+nnB+i]);
/** **/		/*printf("%.2f ",P[nnA+nnB+i]);*/
    }
    fprintf(output,"%.2f %.2f %.2f ", kkoff1, kkoff2, kkcatr);
/** **/ /*printf("%.2f %.2f %.2f ", kkoff1, kkoff2, kkcatr);*/
    fprintf(output,"%.2f %.2f %d %d\n",VV, vv, NOEq, E);
/** **/ /*printf("%.2f %.2f %d %d\n",VV, vv, NOEq, E);*/
	
    fprintf(output,"\n");
    fflush(output);
    fclose(output);
}

void ReadKineticsRates(char *fileName,double *P, double *kkon1, double *kkon2, double *kkcat, double *kkoff1, double *kkoff2, double *kkcatr, int nnA, int nnB, double *VV, double *vv,int *NOEq, int *E)
{
  int i, dummy;
/** **//*printf("read() X\n");*/
  input=fopen(fileName,"r");   
  for(i=0;i<(nnA+nnB);i++)
    dummy=fscanf(input,"%lf ",&P[i]); /*rateT*/
	dummy=fscanf(input, "%lf %lf %lf ", kkon1, kkon2, kkcat);
  for(i=(nnA+nnB+3);i<(nnA+nnB+3+nnA+nnB);i++)
    dummy=fscanf(input,"%lf ",&P[i-3]); /*rateF*/
	dummy=fscanf(input, "%lf %lf %lf\n", kkoff1, kkoff2, kkcatr);
  dummy=fscanf(input,"%lf %lf %d %d\n",VV, vv, NOEq, E);
  fclose(input);
  dummy=dummy;
}

void OLDDomCycleReactIndex(int DC, int *rT, int *rR,int nnA, int nnB)
{
  /** 
   DC: dominant cycle is A(1) or B(2) (DominantCycle)
   *rT: index of toward reaction rate (recT)
   *rR: index of reverse reaction rate (recR)
   nnA: number of the member of cycle A
   nnB: number of the member of cycle B
   **/
    
    if(DC==1)
    {
      *rT=nnA-1;
      *rR=nnA+nnB+nnA-1;
    }
    else if (DC==2)
    {
      *rT=nnA+nnB-1;
      *rR=nnA+nnB+nnA+nnB-1;
    }
    else
    {
      printf("No Dominant cycle -> exit\n");
      exit(1);
    }
}

void DomCycleReactIndex(int DC, int *rT, int *rR,int nnA, int nnB, int PPozRev)
{
  /** 
   DC: dominant cycle is A(1) or B(2) (DominantCycle)
   *rT: index of toward reaction rate (recT)
   *rR: index of reverse reaction rate (recR)
   nnA: number of the member of cycle A
   nnB: number of the member of cycle B
   PPozRev: pozition of revers reactions in the dominant cycle
   **/
/** **//*printf("DomCycleReactIndex\n");*/
    if((DC==1) & (PPozRev>=0))
    {
      *rT=PPozRev/*nnA-1*/;
      *rR=nnA+nnB+PPozRev/*nnA+nnB+nnA-1*/;
    }
    else if ((DC==2) & (PPozRev>=0))
    {
      *rT=nnA+PPozRev/*nnA+nnB-1*/;
      *rR=nnA+nnB+nnA+PPozRev/*nnA+nnB+nnA+nnB-1*/;
    }
    else
    {
      printf("No Dominant cycle -> exit\n");
      exit(1);
    }
    /** **//*printf("X\n");*/
}

