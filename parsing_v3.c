#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#define StrLenPHI 100000

void Write_compactAutcatC(int, double *,double *, double , double, double, char **, int, double*, double, double , double, char *, char **, int*, int, int);
void SimpleRandomNumberInitilaization(gsl_rng *,long int);
void TReaction(char *, char *, char **, int *, int *, int, double*, char *);
void FReaction(char *, char *, char **, int *, int *, int, double*, char *);
void TReaction2(char *, char *, char **, int *, int *, int, double*, char *);
void FReaction2(char *, char *, char **, int *, int *, int, double*, char *);
void WritingGraphViz(char *, char *, char *, int *, int *, double *, double *);

gsl_rng *rcont;

FILE *input, *input2, *output, *fp, *GV;
int main(int argc, char *argv[])
{
	char *file=argv[1]; /** parameter file **/
	int nA=atoi(argv[2]); /** number of member of cycle A **/
	int nB=atoi(argv[3]); /** number of member of cycle B **/
	char line[10000], **eqs, ch[2]="\n" , **eqs2;
	int N, LastTokenLength, NumOfReac, num=0, Tmodel, JoinedCycles, NumOfInteraction, NumOfInteractionT1, NumOfInteractionT2,  kk, ce, *InterActSign=NULL, NOrig, NumOfReacOrig, E, S1, S2, P1, P2, ReacNum=0, komplex1, komplex2, ChemSpec;
  double *InFlowRegime, c, V,v,kon1, koff1, kon2, koff2, kcat, kcatr;	
	char *s=NULL, *ss=NULL,*sss=NULL,*ssss=NULL,PHI[StrLenPHI], sn[1000], FlucEq[100000], **InterAct=NULL;
  int Educt[2]={-1,-1}, Prod[2]={-1,-1},Educt2[3]={-1,-1,-1}, Prod2[3]={-1,-1,-1};
	double *ratesT, *ratesF, *Icc, *Rates, cc, Start, End, DeltaT, MaxCC;
  int LineCounter=0, i=0, j=0, jj=0, index, nindex, jpA, jpB;
  long int InitRandNum=42;
	long long int wcl, OS;
	
  rcont=gsl_rng_alloc(gsl_rng_mt19937); /*gsl_rng_uniform(rcont), gsl_rng_uniform_int(rcont,n)*/
	SimpleRandomNumberInitilaization(rcont,InitRandNum);
  /*strcpy(PHI,"(");*/
	
/** **/	/*printf("file name: %s\n", file);*/

	input=fopen(argv[1],"r");
/***********************************************************/
/** READ General paramaters **/
  if(fgets(line,1000,input)==NULL){;}
/** **/  /*printf("s: %s",line);*/
  ch[0]='\0';
  while(strcmp(ch,"#")!=0)
  {
    if(fscanf(input,"%lf,%lf,%lf,%d,%d,%lf,%lf,%lf,%d,%d,%d,%d\n",&Start, &End, &DeltaT,&NumOfReacOrig,&NOrig, &MaxCC, &V, &v, &Tmodel, &JoinedCycles, &NumOfInteractionT1, &NumOfInteractionT2)!=12){printf("ERROR: fscanf() 1\n");}
/** **/    /*printf("%lf, %lf, %lf, %d, %d, %f, %f, %f, %d, %d, %d, %d\n",Start, End, DeltaT,NumOfReacOrig,NOrig, MaxCC, V, v, Tmodel, JoinedCycles, NumOfInteractionT1, NumOfInteractionT2);*/
    NumOfInteraction = NumOfInteractionT1+NumOfInteractionT2;
    
    ch[0]=fgetc(input);
    fseek(input, ftell(input) - 1, SEEK_SET);
  }
  N=NOrig+(NumOfInteractionT1*1)+(NumOfInteractionT2*2);
	NumOfReac=NumOfReacOrig+(NumOfInteractionT1*2)+(NumOfInteractionT2*3);
/** **//*printf("New NumOfSpecies (N): %d\nNew NumOfReactions (NumOfReac): %d\n", N,NumOfReac);*/
/***********************************************************/
/** **//*exit(120);*/

	ratesT=(double*)calloc(NumOfReac,sizeof(double*));
	ratesF=(double*)calloc(NumOfReac,sizeof(double*));
  Rates=(double*)calloc((2*NumOfReac),sizeof(double*));
	Icc=(double*) calloc(N, sizeof(double));
  InFlowRegime=(double*) calloc(N, sizeof(double));
	for(i=0;i<NumOfReac;i++)
	{
		ratesT[i]=ratesF[i]=-99.99;
/** **/		/*printf("Rates[%d]: %f %f\n",i,ratesT[i],ratesF[i]);	*/	
	}

	
	eqs=(char**)calloc(N,sizeof(char*));
	for(i=0; i<N; i++)
	{
		eqs[i]=(char*)calloc(2*StrLenPHI, sizeof(char));
	}
	
	ss=(char*)calloc(1000, sizeof(char));
	for(i=0; i<N; i++)
	{
		strcpy(eqs[i],"Ith(dydx,");
		sprintf(ss, "%d",i+1);
		strcat(eqs[i], ss);
		strcat(eqs[i], ")=");
/** **/ 		/*printf("%s\n", eqs[i]);*/
    InFlowRegime[i]=999.99;
	}	
	free(ss);

/** **/ /*exit(100);*/
	GV=fopen("GraphRepresentation.dot","w");
	fprintf(GV,"digraph G\n{\n");


  
/** **/ /*fclose(GV);exit(100);*/

/***********************************************************/
/** READ InFlowRegime of certain chemical entities **/
  ch[0]='\0';
  if(fgets(line,1000,input)==NULL){;}
/** **/  /*printf("s: %s",line);*/
  ch[0]=fgetc(input);
  fseek(input, ftell(input) - 1, SEEK_SET);
/** **/  /*printf("ch: %s\n",ch);*/
  index = 0; c=0.0;
  while(strcmp(ch,"#")!=0)
  {
    if(fscanf(input,"%d,%lf\n",&index,&c)!=2){printf("ERROR: fscanf() 2\n");}
/** **/    /*printf("%d,%f\n",index,c);*/
		InFlowRegime[index-1]=c;
/** **/		/*printf("%d,%f\n",index, InFlowRegime[index-1]);*/   
    ch[0]=fgetc(input);
    fseek(input, ftell(input) - 1, SEEK_SET);
  }
  
/** **/  /*for(nindex=1; nindex<=N; nindex++)*/
/** **/		/*printf("%d,%f\n",nindex, InFlowRegime[nindex-1]);*/
  
/** **/ /*exit(1117);*/
/***********************************************************/

/***********************************************************/
/** READ Initial concentrations **/
  ch[0]='\0';
  num = 0;
  if(fgets(line,1000,input)==NULL){;};
/** **/  /*printf("s: %s",line);*/
  ch[0]=fgetc(input);
  fseek(input, ftell(input) - 1, SEEK_SET);
/** ***/  /*printf("ch: %s\n",ch);*/

	while(strcmp(ch,"#")!=0)
  {
    if(fscanf(input,"%d,%lf\n",&index,&cc)!=2){printf("ERROR: fscanf() 3\n");}
/** **/    /*printf("%d,%f\n",index, cc); */  
    if(cc<0.0)
			Icc[index-1]=5.00*gsl_rng_uniform(rcont);/*printf("initial konceentracio (Icc) random feltoltese\n");*/
		else Icc[index-1]=cc;
/** **/ /*printf("Icc[%d]=%f, InFlowRegime[%d]= %f\n", index-1, Icc[index-1], index-1, InFlowRegime[index-1]);*/

    ch[0]=fgetc(input);
    fseek(input, ftell(input) - 1, SEEK_SET);
/** **/ /*printf("ch: %s\n",ch);*/
/** **//*if(index>3 && index<5)exit(115);*/
  }
/** **//*exit(115);*/
/** **/ /*for(j=0;j<N;j++)*/
/** **/ 		/*printf("Icc[%d]=%f, InFlowRegime[%d]= %f\n", j, Icc[j], j, InFlowRegime[j]);*/
/***********************************************************/

/***********************************************************/
/** **//*exit(115);*/
/** READ Reactions **/
	ch[0]='\0';
  num = 0;
  if(fgets(line,1000,input)==NULL){;}
/** **/  /*printf("line: %s",line);*/
  num=0;
	ch[0]=fgetc(input);
  fseek(input, ftell(input) - 1, SEEK_SET);
/** **/  /*printf("ch: %s\n",ch);*/
	
  while(strcmp(ch,"#")!=0)
  {
		s=(char*)calloc(1000, sizeof(char));
		ss=(char*)calloc(1000, sizeof(char));
		sss=(char*)calloc(1000, sizeof(char));
		ssss=(char*)calloc(1000, sizeof(char));

/** **/ /*		printf("%d. sor tokenenként: %s\n",LineCounter, s);*/
		if(fscanf(input,"%d,%d,%lf,%lf,%d,%d\n",&Educt[0], &Educt[1], &ratesT[num], &ratesF[num], &Prod[0], &Prod[1])!=6){printf("ERROR fscanf() 4\n");}
/** **/		/*printf("%d,%d,%f,%f,%d,%d\n",Educt[0], Educt[1], ratesT[num], ratesF[num], Prod[0], Prod[1]);*/
/** **/ /*exit(12);*/
    /*if(ratesT[num]<0)
      ratesT[num]=10.00*gsl_rng_uniform(rcont);
    else;
    if(ratesF[num]<0)
      ratesF[num]=10.00*gsl_rng_uniform(rcont);
    else;*/
    
/** GRAPHVIZ **/    
    WritingGraphViz(ss, sss,ssss, Educt, Prod, &ratesT[0], &ratesF[0]);
		strcpy(sss,"");
		if(Prod[0]==999 && Prod[1]==0) /*Death or Outflow reaction*/
			TReaction(s, ss, eqs, Educt, Prod, num+1, InFlowRegime, sss);
		else
		{
/** **/      /*printf("Normal ág\n");*/
			TReaction(s, ss, eqs, Educt, Prod, num+1, InFlowRegime, sss);
			FReaction(s, ss, eqs, Educt, Prod, num+1, InFlowRegime, sss);
		}
/** **/		/*for(i=0;i<N;i++)*/
/** **/			/*printf("eqs[%d]: %s\n",i,eqs[i]);*/

/** **//*exit(111);*/

		free(s);
		free(ss);
		free(sss);
		free(ssss);

    num++;
		ch[0]=fgetc(input);
    fseek(input, ftell(input) - 1, SEEK_SET);
	}
/** **//*printf("num = %d\n",num);*/
/** **//*exit(17);*/
	ChemSpec=num+2;
	ReacNum=num+1;

/** **/		/*for(i=0;i<N;i++)*/
/** **/			/*printf("eqs[%d]: %s\n",i,eqs[i]);*/
/** **/	/*for(i=0;i<NumOfReac;i++)*/
/** **/    /*printf("Toward[%d]: %f\n",i,ratesT[i]);*/
/** **//*printf("\n\n");*/
/** **//*for(i=0;i<NumOfReac;i++)*/
/** **/   /*printf("Forward[%d]: %f\n",i,ratesF[i]);*/
/** **/ /*exit(111);*/
/***********************************************************/

	
	/*s=(char*)calloc(1000,sizeof(char));
	ss=(char*)calloc(StrLenPHI+100,sizeof(char));

	for(i = 0 ; i < N ; i++)
	{
		strcpy(ss,"-(y");
		sprintf(s,"%d", i+1);
		strcat(ss,s);
		strcat(ss,"*");
    if(MaxCC >0.00)
      strcat(eqs[i],ss);
    else;
	}
	for(i=0;i<N;i++)
			printf("eqs[%d]: %s\n",i,eqs[i]);
	free(s);
	free(ss);*/
/** **/ /*exit(1115);*/ 
/***********************************************************/
/** READ Non-linear interaction Type 1: E+S1=ES1->E+P1+P2 **/
	ch[0]='\0';
	num = 0;
	if(fgets(line,1000,input)==NULL){;}
/** **/ /*printf("line: %s",line);*/
  num=0;
	ch[0]=fgetc(input);
	fseek(input, ftell(input) - 1, SEEK_SET);
/** **/  /*printf("ch: %s\n",ch);*/
/** **/ /*printf("ChemSpec = %d, ReacNum = %d\n",ChemSpec, ReacNum);*/
/** **/ /*exit(17);*/
  while(strcmp(ch,"#")!=0)
	{
    s=(char*)calloc(1000, sizeof(char));
    ss=(char*)calloc(1000, sizeof(char));
		sss=(char*)calloc(1000, sizeof(char));
		ssss=(char*)calloc(1000, sizeof(char));
    
		if(fscanf(input,"%d,%d,%le,%le,%le,%le,%d,%d\n",&S1, &E, &kon1, &koff1, &kcat, &kcatr, &P1,&P2)!=8){printf("ERROR: fscanf() 5\n");}
/** **/			/*printf("%d, %d, %le, %le, %le, %le, %d, %d\n", S1, E, kon1, koff1, kcat, kcatr, P1,P2);*/	
    komplex1 = ChemSpec+(num*1);
/** **/ /*printf("%d + %d = %d -> %d + %d + %d\n", S1, E, komplex1, P1, P2,  E);*/
/** E + S1 = ES1 **/    
    Educt[0] = S1;
    Educt[1] = E;
    Prod[0] = komplex1;
    Prod[1] = 0;
    ratesT[ReacNum+(num*2)-1] = kon1;
    ratesF[ReacNum+(num*2)-1] = koff1;
    TReaction(s, ss, eqs, Educt, Prod, ReacNum+(num*2), InFlowRegime, sss);
    FReaction(s, ss, eqs, Educt, Prod, ReacNum+(num*2), InFlowRegime, sss);

/** **/	/*for(i=0; i<N; i++)*/
/** **/    /*printf("eq[%d]: %s\n",i, eqs[i]);*/
/** **/	/*for(i=0;i<NumOfReac;i++)*/
/** **/    /*printf("Toward[%d]: %f\n",i,ratesT[i]);*/
/** **//*printf("\n\n");*/
/** **//*for(i=0;i<NumOfReac;i++)*/
/** **/   /*printf("Forward[%d]: %f\n",i,ratesF[i]);*/
    
    
/** ES1 - > P1 + P1 +E **/    
    Educt2[0] = komplex1;
    Educt2[1] = 0;
    Educt2[2] = 0;
    Prod2[0] = P1;
    Prod2[1] = P2;
    Prod2[2] = E;
    ratesT[ReacNum+(num*2)] = kcat;
    ratesF[ReacNum+(num*2)] = kcatr;
    TReaction2(s, ss, eqs, Educt2, Prod2, ReacNum+(num*2)+1, InFlowRegime, sss);
    FReaction2(s, ss, eqs, Educt2, Prod2, ReacNum+(num*2)+1, InFlowRegime, sss);
    
/** **/	/*for(i=0; i<N; i++)*/
/** **/    /*printf("%s\n", eqs[i]);*/
/** **/	/*for(i=0;i<NumOfReac;i++)*/
/** **/    /*printf("Toward[%d]: %f\n",i,ratesT[i]);*/
/** **//*printf("\n\n");*/
/** **//*for(i=0;i<NumOfReac;i++)*/
/** **/   /*printf("Forward[%d]: %f\n",i,ratesF[i]);*/

    free(s);
    free(ss);
		free(sss);
		free(ssss);
    
		num++;
		ch[0]=fgetc(input);
		fseek(input, ftell(input) - 1, SEEK_SET);
	}

/** **//*printf("ReacNum= %d, num= %d, komplex1 = %d\n", ReacNum+(num*2), num, komplex1);*/
	if(NumOfInteractionT1>0)
	{
		ChemSpec = NOrig+ (NumOfInteractionT1*1)+1;
		ReacNum= NumOfReacOrig + (NumOfInteractionT1*2)+1;
	}
	else;
	/** **//*printf("ReacNum= %d, num= %d, komplex1 = %d, ChemSpec = %d\n", ReacNum, num, komplex1, ChemSpec);*/
/** **//*exit(16);*/
/***********************************************************/
/** READ Non-linear interaction Type 2: E + S1 = ES1, ES1 + S2 = ES1S2 -> E + P1 + P2 **/
	ch[0]='\0';
	num = 0;
	if(fgets(line,1000,input)==NULL){;}
/** **/ 	/*printf("line: %s",line);*/
  num=0;
	ch[0]=fgetc(input);
	fseek(input, ftell(input) - 1, SEEK_SET);
/** **/  /*printf("ch: %s\n",ch);*/
  while(strcmp(ch,"#")!=0)
	{
    s=(char*)calloc(1000, sizeof(char));
    ss=(char*)calloc(1000, sizeof(char));
		sss=(char*)calloc(1000, sizeof(char));
		ssss=(char*)calloc(1000, sizeof(char));
    
		if(fscanf(input,"%d,%d,%d,%le,%le,%le,%le,%le,%le,%d,%d\n",&S1, &S2, &E, &kon1, &koff1, &kon2, &koff2,&kcat, &kcatr, &P1, &P2)!=11){printf("ERROR: fscanf() 6\n");}
/** **/			/*printf("%d,%d,%d,%le,%le,%le,%le,%le,%le,%d,%d\n",S1, S2, E, kon1, koff1, kon2, koff2,kcat, kcatr, P1, P2);*/	

    komplex1 = ChemSpec+(num*2);
		komplex2 = ChemSpec+(num*2)+1;
/** **/ /*printf("%d + %d = %d, %d + %d = %d -> %d + %d + %d\n", S1, E, komplex1,komplex1, S2,komplex2, P1, P2, E);*/

/** E + S1 = ES1 **/    
    Educt[0] = S1;
    Educt[1] = E;
    Prod[0] = komplex1;
    Prod[1] = 0;
    ratesT[ReacNum+(num*3)-1] = kon1;
    ratesF[ReacNum+(num*3)-1] = koff1;
    TReaction(s, ss, eqs, Educt, Prod, ReacNum+(num*3), InFlowRegime, sss);
    FReaction(s, ss, eqs, Educt, Prod, ReacNum+(num*3), InFlowRegime, sss);

/** **/	/*or(i=0; i<N; i++)*/
/** **/    /*printf("eq[%d]: %s\n",i, eqs[i]);*/
/** **/	/*for(i=0;i<NumOfReac;i++)*/
/** **/    /*printf("Toward[%d]: %f\n",i,ratesT[i]);*/
/** **//*printf("\n\n");*/
/** **//*for(i=0;i<NumOfReac;i++)*/
/** **/   /*printf("Forward[%d]: %f\n",i,ratesF[i]);*/

/** ES1 + S2 = ES1S2 **/    
    Educt[0] = komplex1;
    Educt[1] = S2;
    Prod[0] = komplex2;
    Prod[1] = 0;
    ratesT[ReacNum+(num*3)] = kon2;
    ratesF[ReacNum+(num*3)] = koff2;
    TReaction(s, ss, eqs, Educt, Prod, ReacNum+(num*3)+1, InFlowRegime, sss);
    FReaction(s, ss, eqs, Educt, Prod, ReacNum+(num*3)+1, InFlowRegime, sss);

/** **/	/*for(i=0; i<N; i++)*/
/** **/    /*printf("eq[%d]: %s\n",i, eqs[i]);*/
/** **/	/*for(i=0;i<NumOfReac;i++)*/
/** **/    /*printf("Toward[%d]: %f\n",i,ratesT[i]);*/
/** **//*printf("\n\n");*/
/** **//*for(i=0;i<NumOfReac;i++)*/
/** **/   /*printf("Forward[%d]: %f\n",i,ratesF[i]);*/
    
/** ES1S2 - > P1 + P2 + E **/    
    Educt2[0] = komplex2;
    Educt2[1] = 0;
    Educt2[2] = 0;
    Prod2[0] = P1;
    Prod2[1] = P2;
    Prod2[2] = E;
    ratesT[ReacNum+(num*3)+1] = kcat;
    ratesF[ReacNum+(num*3)+1] = kcatr;
    TReaction2(s, ss, eqs, Educt2, Prod2, ReacNum+(num*3)+2, InFlowRegime, sss);
    FReaction2(s, ss, eqs, Educt2, Prod2, ReacNum+(num*3)+2, InFlowRegime, sss);
    
/** **/	/*for(i=0; i<N; i++)*/
/** **/    /*printf("%s\n", eqs[i]);*/
/** **/	/*for(i=0;i<NumOfReac;i++)*/
/** **/    /*printf("Toward[%d]: %f\n",i,ratesT[i]);*/
/** **//*printf("\n\n");*/
/** **//*for(i=0;i<NumOfReac;i++)*/
/** **/   /*printf("Forward[%d]: %f\n",i,ratesF[i]);*/
/** **//*exit(15);*/
    free(s);
    free(ss);
		free(sss);
		free(ssss);
    
		num++;
		ch[0]=fgetc(input);
		fseek(input, ftell(input) - 1, SEEK_SET);
	}	
	
/** **//*exit(133);*/


/***********************************************************/
/** READ JoinedCycles at the member of cycles **/
	ch[0]='\0';
	num = 0;
	if(fgets(line,1000,input)==NULL){;}
/** */	/*printf("line: %s",line);*/
	num=0;
	ch[0]=fgetc(input);
	fseek(input, ftell(input) - 1, SEEK_SET);
/** **/  /*printf("ch: %s\n",ch);*/
	while(feof(input)==0)
	{
		if(fscanf(input,"%d=%d\n",&jpA, &jpB)!=2){printf("ERROR: fscanf() 7\n");}
		num++;
		/*ch[0]=fgetc(input);
		fseek(input, ftell(input) - 1, SEEK_SET);*/
	}
	/*if(JoinedCycles==0)
		FlucEq[0]='\0';*/
/** **//*printf("jpA = %d, jpB = %d\n",jpA, jpB);*/
/***********************************************************/
/** **//*xit(133);*/

	
/** **//*exit(1230);*/

/***********************************************************/

  fclose(input);
  
/** **/ /*exit(1115); */
	s=(char*)calloc(1000,sizeof(char));
/*
 If no cemostat, the X regime may be 
 *1) unlimited: dX/dt=c where c=0; 
 *2) constant flow dX/dt= c, where c>0
 *3) limited dx/dt=c, where c= -k*A*B + ....
*/

	if(Tmodel==0)
	{
		for(i=0; i<N; i++)
		{
			if(InFlowRegime[i]>=0.00 && InFlowRegime[i]<999.99)
			{
				strcpy(eqs[i],"Ith(dydx,");
				sprintf(s, "%d",i+1);
				strcat(eqs[i], s);
				strcat(eqs[i], ")=");
				sprintf(s,"%f",InFlowRegime[i]);
				strcat(eqs[i], s);
			}
		}
	}
	
	free(s);

	s=(char*)calloc(1000,sizeof(char));

	if(jpB>0)
				strcat(eqs[jpB-1], "0.00");
	if(Tmodel>0)
	{
		for(i=0; i<N; i++)
		{
			if((i+1)==jpB)
				strcat(eqs[i], " ");
			else
			{
				strcat(eqs[i], "+((v/V)*(cy");
				sprintf(s,"%d", i+1);
				strcat(eqs[i],s);
				strcat(eqs[i], "-y");
				strcat(eqs[i],s);
				strcat(eqs[i],"))");
			}
		}
	}
	
	free(s);
	
	for(i=0; i<N; i++)
		strcat(eqs[i], ";");
	
/** **/	/*for(i=0; i<N; i++)*/
/** **/    /*printf("%s\n", eqs[i]);*/
	
/** **/	/*for(j=0;j<N;j++)*/
/** **/		/*printf("Icc[%d]=%f, InFlowRegime[%d]= %f\n", j, Icc[j], j, InFlowRegime[j]);*/
/** **//*printf("\n");*/
	
	for(i=0;i<NumOfReac;i++)
  {
    Rates[i]=ratesT[i];
		Rates[i+NumOfReac]=ratesF[i];
/** **/     /*printf("Toward[%d]: %f\tForward[%d]: %f\n",i,Rates[i],i+NumOfReac,Rates[i+NumOfReac]);*/
  }
	
	for(i=0 ; i< nA ; i++)
	{
		if((i+1)==jpA)
			fprintf(GV, "%d [shape=box, color=green, style=bold, label=\"A%d=B%d\"];\n", i+1, i+1, jpB-nA);
		else fprintf(GV, "%d [shape=box, color=red, style=bold, label=A%d];\n", i+1, i+1);
	}
	if(nB>0)
	{
		for(i=nA ; i< nA+nB ; i++)
		{
			if((i+1)==jpB)
				fprintf(GV, "/*%d [shape=box, color=blue, style=bold, label=B%d];*/\n", i+1, (i+1)-nA);
			else fprintf(GV, "%d [shape=box, color=blue, style=bold, label=B%d];\n", i+1, (i+1)-nA);
		}
		fprintf(GV, "%d [shape=box, color=black, style=bold, label=X];\n", nA+nB+1);

	}
	else
		fprintf(GV, "%d [shape=box, color=black, style=bold, label=X];\n", nA+1);
	fprintf(GV,"}\n");
	fclose(GV);
  
/** **/ /*exit(100);*/
  
  Write_compactAutcatC(N, &Icc[0], &Rates[0], Start, End, DeltaT, eqs, NumOfReac, &InFlowRegime[0], V, v, MaxCC, FlucEq, &InterAct[0], &InterActSign[0], NumOfInteraction, jpB);

   for(i=0; i<N; i++)
    free(eqs[i]);
   free(eqs);
   free(Rates);
   free(Icc);
   free(InFlowRegime);
   free(ratesF);
   free(ratesT);
   gsl_rng_free (rcont);
   
   return 0;
}
/** C kódok generálása **/	
void Write_compactAutcatC(int NN,  double *IIcc, double *PPars, double SStart, double EEnd, double DDeltaT, char **Eqs, int NNumOfReac, double *IInFlowRegime, double VV, double vv, double MMaxCC, char *FFlucEq, char ** IInterAct, int *IInteractSign, int NNumOfInteraction, int jjpB)
{
  int i;

 /** cv_ODE_solver.c **/ 
  output=fopen("compactAutCat.c","wt");
  fprintf(output,"#include <stdio.h>\n");
  fprintf(output,"#include <stdlib.h>\n");
  fprintf(output,"#include <math.h>\n");
  fprintf(output,"#include <time.h>\n");
  fprintf(output,"#include <string.h>\n\n");
  fprintf(output,"#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */\n");
  fprintf(output,"#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */\n");
  fprintf(output,"#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix    */\n");
  fprintf(output,"#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver*/\n");
  fprintf(output,"#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype*/\n\n\n");
  fprintf(output,"/** **/ /* absolute error */ \n");
  fprintf(output,"#define ABSTOL RCONST(1.0e-10)\n");
  fprintf(output,"/** **/ /* relative error */ \n");
  fprintf(output,"#define RELTOL RCONST(1.0e-8) \n");
  fprintf(output,"#define NOUT 20000000000\n");
  fprintf(output,"#define Ith(v,i) NV_Ith_S(v,i-1)       /* Ith numbers components 1..Numbers Of Equations  */\n");
	fprintf(output,"#define TOL (1.0e-6)\n\n\n");
  fflush(output);

	fprintf(output,"int NumOfVars=%d, NumOfPars=%d;\n",NN, (2*NNumOfReac)+2);
	fprintf(output,"double InitConc[%d]={",NN);
	for(i=0;i<NN-1; i++)
    fprintf(output,"%f,",IIcc[i]);
  fprintf(output,"%f};\n", IIcc[NN-1]);	
  fprintf(output,"double Pars[%d]={",(2*NNumOfReac)+2);
	for(i=0;i<(2*NNumOfReac); i++)
    fprintf(output,"%f,",PPars[i]);
	fprintf(output,"%f,",VV);
	fprintf(output,"%f",vv);
  fprintf(output,"};\n");
/** **//*printf("1X\n"); */ 
/** **//*printf("PPars[%d]=%f\n",2*NNumOfReac,PPars[(2*NNumOfReac)]);*/
  /*fprintf(output,"%f};\n", PPars[(2*NNumOfReac)]);	
printf("2X\n");*/
  fprintf(output,"realtype START=%f, END = %f, StepT=%f;\n  ",SStart, EEnd, DDeltaT);
  fprintf(output,"FILE  *output,*input, *fp, *output2;\n\n");
  fflush(output);
/** **//*printf("3X\n");*/
  fprintf(output,"int eqs(realtype, N_Vector , N_Vector, void *);\n");
  fprintf(output,"static int check_flag(void *flagvalue, char *funcname, int opt);\n\n");

/** FUNCTION main() **/
	fprintf(output,"int main(int argc, char *argv[])\n");
  fprintf(output,"{\n\n");
  
	fprintf(output,"	int nA=atoi(argv[1]);\n");
	fprintf(output,"	int nB=atoi(argv[2]);\n");
  fprintf(output,"	int j,flag, iout, num, /*num2, iout2, num3,*/ type,jpB=%d;\n", jjpB);
  fprintf(output,"	double InitSumCC = 0.00/*, SumCC = 0.00*/, *OriginalPars=NULL, CC[nA+nB][3], SumA, SumB");
  
  /*if(NNumOfInteraction>0)
	{
		fprintf(output,",D[%d]={",NNumOfInteraction);
		for(i=0;i<NNumOfInteraction;i++)
		{
			if(i==NNumOfInteraction-1)
				fprintf(output,"atof(argv[%d])}",i+1);
			else fprintf(output,"atof(argv[%d]),",i+1);
		}
		fprintf(output,", K[%d]={",NNumOfInteraction);
		for(i=0;i<NNumOfInteraction;i++)
		{
			if(i==NNumOfInteraction-1)
				fprintf(output,"atof(argv[%d])}",i+NNumOfInteraction+1);
			else fprintf(output,"atof(argv[%d]),",i+NNumOfInteraction+1);
		}
	}
	else;*/
	fprintf(output,";\n");
  fprintf(output,"	realtype reltol = RELTOL,  T, nextT;\n");
  fprintf(output,"	N_Vector Conc = NULL, abstol = NULL;\n");
  fprintf(output,"	SUNMatrix A;\n");
  fprintf(output,"	SUNLinearSolver LS;\n");
  fprintf(output,"	void *cvode_mem = NULL;\n\n");
	
	
	fprintf(output,"	SumA=SumB=0.00;\n");
	fprintf(output,"	for(j=0;j<nA;j++)\n");
	fprintf(output,"		SumA+=InitConc[j];\n");
	fprintf(output,"	for(j=nA;j<(nA+nB);j++)\n");
	fprintf(output,"		SumB+=InitConc[j];\n");
	fprintf(output,"	for(j=0;j<(nA+nB);j++)\n");
	fprintf(output,"	{\n");
	fprintf(output,"		CC[j][0]=InitConc[j];\n");
	fprintf(output,"		CC[j][1]=0.00;\n");
	fprintf(output,"		CC[j][2]=0.00;\n");
	fprintf(output,"	}\n\n");
	
	
  if(MMaxCC>0.00)
  {
    /*fprintf(output,"	double ");
    for(i=0;i<NNumOfReac;i++)
      fprintf(output, "  k%d = Pars[%d], ",i+1,i);
    for(i=NNumOfReac;i<(2*NNumOfReac)-1;i++)
      fprintf(output, "  k%dr = Pars[%d], ",i-NNumOfReac+1,i);
    fprintf(output, "  k%dr = Pars[%d];\n\n",(2*NNumOfReac)-NNumOfReac,(2*NNumOfReac)-1);*/
  }
  
  fprintf(output,"	OriginalPars=(double*)calloc(NumOfPars, sizeof(double));\n");
	/*fprintf(output,"	output2=fopen(\"k_temp\",\"w\");\n");*/
  fprintf(output,"	for(j=0;j<NumOfPars;j++)\n");
	fprintf(output,"	{\n");
  fprintf(output,"		OriginalPars[j]=Pars[j];\n");
	/*fprintf(output,"		fprintf(output2,\"%%f \",Pars[j]);\n");*/
	fprintf(output,"	}\n");
	/*fprintf(output,"	fprintf(output2,\"\\n\");\n");*/
	/*fprintf(output,"	fclose(output2);\n");*/
	fprintf(output,"	Conc = N_VNew_Serial(NumOfVars);\n");
  fprintf(output,"	if (check_flag((void *)Conc, \"N_VNew_Serial\", 0)) return(1);\n");
  fprintf(output,"	abstol = N_VNew_Serial(NumOfVars);\n");
  fprintf(output,"	if (check_flag((void *)abstol, \"N_VNew_Serial\", 0)) return(1);\n\n");
  fflush(output);  
  
/** **/	//fprintf(output,"/** **/  fr=(int)((END/StepT)/%f);/*mikor kell váltani a v-ben*/\n", FlucFreq);*/
  fprintf(output,"	InitSumCC = 0.00;\n");
  fprintf(output,"	output=fopen(\"res.dat\", \"w\");/* external file name */\n");
  //if(MMaxCC>0.00)
  //  fprintf(output,"  fp=fopen(\"res2.dat\", \"w\");/* external file name */\n");
  fprintf(output,"	fprintf(output,\"%%e \",START);\n");
  //if(MMaxCC>0.00)
    //fprintf(output,"  fprintf(fp,\"%%e \",START);\n");
  fprintf(output,"	for(j=0;j<NumOfVars;j++)\n");
  fprintf(output,"	{\n"); 
  fprintf(output,"		fprintf(output,\"%%e \",InitConc[j]);\n");
  //if(MMaxCC>0.00)
    //fprintf(output,"    fprintf(fp,\"%%e \",0.000000);\n");
  fprintf(output,"		InitSumCC+=InitConc[j];\n");
  fprintf(output,"		Ith(Conc,j+1) = InitConc[j];\n");
  fprintf(output,"		Ith(abstol,j+1) = ABSTOL;\n");
  fprintf(output,"	}\n");
  fprintf(output,"	fprintf(output, \"\\n\");\n");
  fprintf(output,"	fclose(output);\n");
	fprintf(output,"	fflush(output);\n");
  //if(MMaxCC>0.00)
  //{
   // fprintf(output,"  fprintf(fp, \"\\n\");\n");
    //fprintf(output,"  fclose(fp);\n\n\n");
  //}	
  fprintf(output,"	\n\n");
	
  
  fprintf(output,"	cvode_mem = CVodeCreate(CV_BDF); /* CV_BDF: stiff problems, CV_ADAMS: non stiff problems */\n");
  fprintf(output,"	if (check_flag((void *)cvode_mem, \"CVodeCreate\", 0)) return(1);\n");
  fprintf(output,"	flag = CVodeSetUserData(cvode_mem, &Pars[0]);\n");
  fprintf(output,"	if (check_flag(&flag, \"CVodeSetUserData\", 1)) return(1);\n");
  fprintf(output,"	flag = CVodeInit(cvode_mem, eqs, START, Conc);\n");
  fprintf(output,"	if (check_flag(&flag, \"CVodeInit\", 1)) return(1);\n");
  fprintf(output,"	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);\n");
  fprintf(output,"	if (check_flag(&flag, \"CVodeSVtolerances\", 1)) return(1);\n");
  fprintf(output,"	A = SUNDenseMatrix(NumOfVars, NumOfVars);/* Create dense SUNMatrix for use in linear solves */\n");
  fprintf(output,"	if(check_flag((void *)A, \"SUNDenseMatrix\", 0)) return(1);\n");
  fprintf(output,"	LS = SUNLinSol_Dense(Conc, A);/* Create dense SUNLinearSolver object for use by CVode */\n");
  fprintf(output,"	if(check_flag((void *)LS, \"SUNLinSol_Dense\", 0)) return(1);\n");
  fprintf(output,"	flag = CVodeSetLinearSolver(cvode_mem, LS, A);/* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */\n\n\n");
  fflush(output);
    
  fprintf(output,"	nextT = START+StepT;\n");
  fprintf(output,"	iout = 0;\n");
	fprintf(output,"	/*iout2=-99;*/\n");
	fprintf(output,"	/*num3=num2=0;*/\n");
  fprintf(output,"	while(1)\n");
  fprintf(output,"	{\n");
  fprintf(output,"		flag = CVode(cvode_mem, nextT, Conc, &T, CV_NORMAL);\n");
  fprintf(output,"		/*num=0;*/\n");
	fprintf(output,"		output=fopen(\"res.dat\", \"a\");/* external file name */\n");
	fprintf(output,"/** **/    /*printf(\"%%e \",T);*/\n");
	fprintf(output,"		fprintf(output,\"%%e \",nextT);\n");
  fprintf(output,"		for(j=0;j<NumOfVars;j++)\n");
  fprintf(output,"		{\n");
	fprintf(output,"/** **/    	 /*printf(\"%%e \",Ith(Conc,j+1));*/\n");
  fprintf(output,"		fprintf(output,\"%%e \",Ith(Conc,j+1));\n");
  fprintf(output,"		}\n");
	fprintf(output,"/** **//*printf(\"\\n\");*/\n");
	fprintf(output,"		fprintf(output,\"\\n\");\n");
	fprintf(output,"		fflush(output);\n");
	fprintf(output,"		fclose(output);\n");
	
	fprintf(output,"		type = 0 ;\n");
	fprintf(output,"		SumA=SumB=0.0;\n");
  fprintf(output,"		if(iout%1000==0)\n");
  fprintf(output,"		{\n");
	fprintf(output,"			for(j = 0; j< nA ; j++)\n");
	fprintf(output,"				SumA+=Ith(Conc,j+1);\n");
	fprintf(output,"			if(nB>0)\n");
	fprintf(output,"			{\n");
	fprintf(output,"				for(j = nA; j< (nA+nB) ; j++)\n");
	fprintf(output,"					SumB+=Ith(Conc,j+1);\n\n");
		
	fprintf(output,"/** **//*printf(\"t= %%d: num2= %%d, num3=%%d, \", iout, num2, num3);*/\n");			
	fprintf(output,"				if(SumA<TOL && SumB<TOL)\n");
	fprintf(output,"					type=8; /** Cycle A and B died out **/\n");
	fprintf(output,"				if(SumA<TOL && SumB>TOL)\n");
	fprintf(output,"					type=2; /** Winner : B **/\n");
	fprintf(output,"				if(SumA>TOL && SumB<TOL)\n");
	fprintf(output,"					type=1; /** Winner : A **/\n");
	fprintf(output,"			}\n");
	fprintf(output,"			else\n");
	fprintf(output,"			{\n");
	fprintf(output,"				if(SumA<TOL )\n");
	fprintf(output,"					type=8; /** Cycle A died out **/\n");
	fprintf(output,"			}\n");
	
	fprintf(output,"			/*else if(num2>0 && num3>0)\n");
	fprintf(output,"			type=4;*/ /** System died **/\n");
	fprintf(output,"			/*else;*//* type = 9;*/ /** PROBLEM **/\n\n");
	fprintf(output,"			num=0;\n");
	fprintf(output,"			if(type==0)\n");
	fprintf(output,"			{\n");
	fprintf(output,"				for(j = 0; j< (nA+nB) ; j++)\n");
	fprintf(output,"				{\n");
	fprintf(output,"					if((j+1)==jpB)\n");
	fprintf(output,"						num++; \n");
	fprintf(output,"					else\n");
	fprintf(output,"					{\n");
	fprintf(output,"						CC[j][1]=Ith(Conc,j+1);\n");
	fprintf(output,"						CC[j][2]=fabs((CC[j][0]-CC[j][1])/CC[j][1]);\n");
	fprintf(output,"						if(CC[j][2]<TOL)\n");
	fprintf(output,"							num++;\n");
	fprintf(output,"						else;\n");
	fprintf(output,"					}\n");
	fprintf(output,"				}\n");
		
	fprintf(output,"				if(num==(nA+nB))\n");
	fprintf(output,"					type = 3;\n\n");
	fprintf(output,"				for(j = 0; j< (nA+nB) ; j++)\n");
	fprintf(output,"				{\n");
	fprintf(output,"					CC[j][0]=CC[j][1];\n");
	fprintf(output,"					CC[j][2]=0.00;\n");
	fprintf(output,"				}\n");
	fprintf(output,"/** **/ /*printf(\"%%d\\n\", num2);*/\n\n");
	fprintf(output,"			}\n");

	fprintf(output,"			if(type>0)\n");
	fprintf(output,"				return type;\n");
  fprintf(output," 		}\n");
	
	//if(MMaxCC>0.00)
  //{
   // fprintf(output,"    fp=fopen(\"res2.dat\", \"a\");/* external file name */\n");
    //fprintf(output,"    fprintf(fp,\"%%e \",nextT);\n");
    //for(i=0;i<NN; i++)    
    //{
     // fprintf(output,"    fprintf(fp,\"%%e \",");
      //fprintf(output,");\n");
    //}
    //fprintf(output,"    fprintf(fp,\"\\n\");\n");
    //fprintf(output,"    fclose(fp);\n\n");
  //}
  fflush(output);  
  
	if(strlen(FFlucEq)!=0)
	{
		fprintf(output,"		num=1;\n");
		fprintf(output,"    if(num>0)\n");
		fprintf(output,"		{\n");
		fprintf(output,"			Pars[%d]= %s;\n",(2*NNumOfReac)+1, FFlucEq);
		fprintf(output,"      flag = CVodeReInit(cvode_mem, nextT, Conc);\n\n");
		fprintf(output,"			/*printf(\"T: %%f, Valtozaskor: v= %%f\\n\", T,Pars[%d]);*/\n",(2*NNumOfReac)+1);
	fprintf(output,"						}\n");
	}
	/*else if(NNumOfInteraction>0)
	{
		fprintf(output,"	num=1;\n");
		fprintf(output,"	if(num>0)\n");
		fprintf(output,"	{\n");
		for(i=0;i<NNumOfInteraction;i++)
			fprintf(output,"		%s;\n",IInterAct[i]);
		fprintf(output,"		flag = CVodeReInit(cvode_mem, nextT, Conc);\n\n");
	fprintf(output,"	}\n");
	}*/
	else
	{
		;/*fprintf(output,"    if(num>0)\n");
		fprintf(output,"      flag = CVodeReInit(cvode_mem, nextT, Conc);\n\n");*/
	}
	fprintf(output,"\n");
  
	fprintf(output,"		if (check_flag(&flag, \"CVode\", 1)) break;\n");
  fprintf(output,"		if (flag == CV_SUCCESS)\n");
  fprintf(output,"		{\n");
  fprintf(output,"			iout++;\n");
  fprintf(output,"			nextT += StepT/*StepT*/;\n");
  fprintf(output,"		}\n");
	fprintf(output,"/** **//*printf(\"iout= %%d, iout2= %%d, num3 = %%d\\n\",iout, iout2, num3);	*/\n");
  fprintf(output,"		if (/*iout == NOUT ||*/ T >= END /*|| iout==iout2*/) break; /*running time: NOUT*StepT*/\n");
  fprintf(output,"	}\n\n");
  
  fprintf(output,"	N_VDestroy(Conc); /* Free Conc vector */\n");
  fprintf(output,"	N_VDestroy(abstol); /* Free abstol vector */  \n");
  fprintf(output,"	CVodeFree(&cvode_mem); /* Free integrator memory */\n");
  fprintf(output,"	SUNLinSolFree(LS); /* Free the linear solver memory */\n");
  fprintf(output,"	SUNMatDestroy(A); /* Free the matrix memory */\n\n");
  
	fprintf(output,"	free(OriginalPars); \n\n");
	
	fprintf(output,"  type=9;\n");
  fprintf(output,"  return type;\n");
  fprintf(output,"}\n\n\n");
  fflush(output);
  
/** FUNCTION check_flag() **/
  fprintf(output,"static int check_flag(void *flagvalue, char *funcname, int opt)\n");
  fprintf(output,"{\n");
  fprintf(output,"  int *errflag;\n\n");
  
  fprintf(output,"  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */\n");
  fprintf(output,"  if (opt == 0 && flagvalue == NULL) {\n");
  fprintf(output,"    fprintf(stderr, \"\\nSUNDIALS_ERROR: %%s() failed - returned NULL pointer\\n\\n\",funcname);\n");
  fprintf(output,"    return(1);\n");
  fprintf(output,"  }\n\n");
  
  fprintf(output,"  /* Check if flag < 0 */\n");
  fprintf(output,"  else if (opt == 1) {\n");
  fprintf(output,"    errflag = (int *) flagvalue;\n");
  fprintf(output,"    if (*errflag < 0) {\n");
  fprintf(output,"       fprintf(stderr, \"\\nSUNDIALS_ERROR: %%s() failed with flag = %%d\\n\\n\",funcname, *errflag);\n");
  fprintf(output,"  return(1); }}\n\n");
  
  fprintf(output,"  /* Check if function returned NULL pointer - no memory allocated */\n");
  fprintf(output,"  else if (opt == 2 && flagvalue == NULL) {\n");
  fprintf(output,"    fprintf(stderr, \"\\nMEMORY_ERROR: %%s() failed - returned NULL pointer\\n\\n\",funcname);\n");
  fprintf(output,"    return(1); }\n\n");
  
  fprintf(output,"  return(0);\n");
  fprintf(output,"}\n\n\n");
  fflush(output);
  
/** FUNCTION eqs() **/
  fprintf(output,"int eqs(realtype t, N_Vector Conc, N_Vector dydx, void *user_data)\n");
  fprintf(output, "{\n");
  fprintf(output,"  realtype ");
  for(i=0;i<NN;i++)
    fprintf(output, "y%d, ", i+1);

  for(i=0;i<NNumOfReac;i++)
    fprintf(output, " k%d, ", i+1);

  for(i=0;i<NNumOfReac;i++)
    fprintf(output, " k%dr, ", i+1);

  if(VV>0.00)
  {
	for(i=0;i<NN;i++)
		fprintf(output, " cy%d, ", i+1);
  }
  fprintf(output, "V,v;\n");
  fprintf(output, "  double *p;\n");
  fprintf(output, "  p=(double*)user_data;\n");
/** értékadás **/	
  for(i=0;i<NN;i++)
    fprintf(output, "  y%d = Ith(Conc, %d);\n",i+1,i+1);
  for(i=0;i<NNumOfReac;i++)
    fprintf(output, "  k%d = p[%d];\n",i+1,i);
  for(i=NNumOfReac;i<(2*NNumOfReac);i++)
    fprintf(output, "  k%dr = p[%d];\n",i-NNumOfReac+1,i);

	
  if(VV>0.00)
  {
	fprintf(output,"	V = p[%d];\n",2*NNumOfReac);
	fprintf(output,"	v = p[%d];\n",(2*NNumOfReac)+1);
	for(i=0;i<NN;i++)
	{
		if(IInFlowRegime[i]==999.99)
		fprintf(output, "  cy%d = %f;\n",i+1,0.00);
		else fprintf(output, "  cy%d = %f;\n",i+1,IInFlowRegime[i]);
	}
  }
  /*fprintf(output, "  c = p[%d];\n\n",2*NNumOfReac);*/
  fflush(output);
 
	fprintf(output,"	/*printf(\"Integralaskor: v= %%f(%%f)\\n\", p[%d],v);*/\n",(2*NNumOfReac)+1);
  for(i=0;i<NN;i++)
      fprintf(output, "  %s\n",Eqs[i]);
  fflush(output);
 
  fprintf(output,"\n  return(0);\n");
  fprintf(output,"}\n\n");
  fflush(output);
  fclose(output);
  
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
		
		gsl_rng_set(r,InitVOfRN);
	}
	/** **/ /*printf("InitVOfRN= %ld, InitialRandNumber: %ld\n",InitVOfRN,RandNum);*/
	/** **/ /*printf("### END SimpleRandomNumberInitilaization()\n");*/
}


void TReaction(char *S, char *SS, char **EQS, int *Ed, int *Pr, int ind, double *InFlowREGIME, char *SSS)
{
	int j;
/** **//*printf("TReaction: %d\n", ind);*/
	/** TOWARD REACTION**/			
		strcpy(S,"(k");
		sprintf(SS,"%d",ind);
		strcat(S,SS);
		strcat(S,"*");
			
		if(Ed[0]!=0)
		{
			strcat(S,"y");
			sprintf(SS,"%d",Ed[0]);
			strcat(S,SS);	
		}
		else 
		{
			strcat(S,"");
		}
/** **/    /*printf("ED0: %s\n",S);*/
			
		if(Ed[1]!=0)
		{
			strcat(S,"*y");
			sprintf(SS,"%d",Ed[1]);
			strcat(S,SS);	
			strcat(S,")");
		}
		else 
		{
			strcat(S,")");
		}
/** **/ 		/*printf("ED1: %s\n",S);*/
    
		for(j = 0; j<2;j++)
		{
/** **/      /*printf("Educt(%d): %d: %f\n",j,Ed[j],InFlowREGIME[Ed[j]-1]);*/
/** **/      /*printf("Educt(%d): %d\n",j,Ed[j]);*/
			if(Ed[j]!=0)
			{				
					strcat(EQS[Ed[j]-1],"-");
					strcat(EQS[Ed[j]-1],S);
			}
/** **/      /*printf("Prod(%d): %d: %f\n",j,Pr[j],InFlowREGIME[Pr[j]-1]);*/
/** **/      /*printf("Prod(%d): %d\n",j,Pr[j]);*/
			if(Pr[j]!=0 && Pr[j]!=999)
			{
				strcat(EQS[Pr[j]-1],"+");
				strcat(EQS[Pr[j]-1],S);
			}
		}
		

}

void FReaction(char *S, char *SS, char **EQS, int *Ed, int *Pr, int ind, double *InFlowREGIME, char *SSS)
{
	int j;
/** FORWARD REACTION**/	
/** **//*printf("FReaction: %d\n", ind);*/
/** **//*printf("%s\n",S);*/
		strcpy(S,"(k");
/** **/ /*printf("%s\n",S);*/
		sprintf(SS,"%d",ind);
		strcat(S,SS);
		strcat(S,"r*");
			
		if(Pr[0]!=0)
		{
			strcat(S,"y");
			sprintf(SS,"%d",Pr[0]);
			strcat(S,SS);			
		}
		else 
		{
			strcat(S,"");
		}
/** **/    /*printf("Pr0: %s\n",S);*/

		if(Pr[1]!=0)
		{
			strcat(S,"*y");
			sprintf(SS,"%d",Pr[1]);
			strcat(S,SS);	
			strcat(S,")");
		}
		else 
		{
			strcat(S,")");
		}
/** **/    /*printf("Pr1: %s\n",S);*/

		for(j = 0; j<2;j++)
		{
/** **/			/*printf("Educt: %d: %f\n",Ed[j],InFlowREGIME[Ed[j]-1]);*/
/** **/      /*printf("Educt(%d): %d\n",j,Ed[j]);*/
			if(Ed[j]!=0)
			{
					strcat(EQS[Ed[j]-1],"+");
					strcat(EQS[Ed[j]-1],S);
			}
/** **/			/*printf("Prod(%d): %d\n",j,Pr[j]);*/
			if(Pr[j]!=0)
			{
				strcat(EQS[Pr[j]-1],"-");
				strcat(EQS[Pr[j]-1],S);
			}
		}
}

void TReaction2(char *S, char *SS, char **EQS, int *Ed, int *Pr, int ind, double *InFlowREGIME, char *SSS)
{
	int j;
/** **//*printf("TReaction: %d\n", ind);*/
	/** TOWARD REACTION**/			
		strcpy(S,"(k");
		sprintf(SS,"%d",ind);
		strcat(S,SS);
		strcat(S,"*");
			
		if(Ed[0]!=0)
		{
			strcat(S,"y");
			sprintf(SS,"%d",Ed[0]);
			strcat(S,SS);	
		}
		else 
		{
			strcat(S,"");
		}
/** **/    /*printf("ED0: %s\n",S);*/
			
		if(Ed[1]!=0)
		{
			strcat(S,"*y");
			sprintf(SS,"%d",Ed[1]);
			strcat(S,SS);	
		}
		else 
		{
			strcat(S,"");
		}
/** **/    /*printf("ED1: %s\n",S);*/	
			
		if(Ed[2]!=0)
		{
			strcat(S,"*y");
			sprintf(SS,"%d",Ed[2]);
			strcat(S,SS);	
			strcat(S,")");
		}
		else 
		{
			strcat(S,")");
		}
/** **/ 		/*printf("ED2: %s\n",S);*/
    
		for(j = 0; j<3;j++)
		{
/** **/      /*printf("Educt(%d): %d: %f\n",j,Ed[j],InFlowREGIME[Ed[j]-1]);*/
/** **/      /*printf("Educt(%d): %d\n",j,Ed[j]);*/
			if(Ed[j]!=0)
			{				
					strcat(EQS[Ed[j]-1],"-");
					strcat(EQS[Ed[j]-1],S);
			}
/** **/      /*printf("Prod(%d): %d: %f\n",j,Pr[j],InFlowREGIME[Pr[j]-1]);*/
/** **/      /*printf("Prod(%d): %d\n",j,Pr[j]);*/
			if(Pr[j]!=0 && Pr[j]!=999)
			{
				strcat(EQS[Pr[j]-1],"+");
				strcat(EQS[Pr[j]-1],S);
			}
		}
		

}

void FReaction2(char *S, char *SS, char **EQS, int *Ed, int *Pr, int ind, double *InFlowREGIME, char *SSS)
{
	int j;
/** FORWARD REACTION**/	
/** **//*printf("FReaction: %d\n", ind);*/
/** **//*printf("%s\n",S);*/
		strcpy(S,"(k");
/** **/ /*printf("%s\n",S);*/
		sprintf(SS,"%d",ind);
		strcat(S,SS);
		strcat(S,"r*");
			
		if(Pr[0]!=0)
		{
			strcat(S,"y");
			sprintf(SS,"%d",Pr[0]);
			strcat(S,SS);			
		}
		else 
		{
			strcat(S,"");
		}
/** **/    /*printf("Pr0: %s\n",S);*/

    if(Pr[1]!=0)
		{
			strcat(S,"*y");
			sprintf(SS,"%d",Pr[1]);
			strcat(S,SS);			
		}
		else 
		{
			strcat(S,"");
		}
/** **/    /*printf("Pr1: %s\n",S);*/
		
		if(Pr[2]!=0)
		{
			strcat(S,"*y");
			sprintf(SS,"%d",Pr[2]);
			strcat(S,SS);	
			strcat(S,")");
		}
		else 
		{
			strcat(S,")");
		}
		
/** **/    /*printf("Pr2: %s\n",S);*/

		for(j = 0; j<3;j++)
		{
/** **/			/*printf("Educt: %d: %f\n",Ed[j],InFlowREGIME[Ed[j]-1]);*/
/** **/      /*printf("Educt(%d): %d\n",j,Ed[j]);*/
			if(Ed[j]!=0)
			{
					strcat(EQS[Ed[j]-1],"+");
					strcat(EQS[Ed[j]-1],S);
			}
/** **/			/*printf("Prod(%d): %d\n",j,Pr[j]);*/
			if(Pr[j]!=0)
			{
				strcat(EQS[Pr[j]-1],"-");
				strcat(EQS[Pr[j]-1],S);
			}
		}
}



void WritingGraphViz(char *SS, char *SSS, char *SSSS, int *Ed, int *Pr, double *RatesT, double *RatesF)
{
    int j,jj; 
    char sp[100];
/** **/    /*printf("#### In WritingGraphViz() ####\n");*/
    for(j=0;j<2;j++)
    {      
      strcpy(SSS,"");
      if(Ed[j]!=0)
      {
        sprintf(SS,"%d",Ed[j]);
        strcat(SSS,SS);
        strcat(SSS,"->");
/** **/        /*printf("GV1: %s\n",SSS);*/
        for(jj=0;jj<2;jj++)
        {
					strcpy(SSSS,"");
          if(Pr[jj]!=0 && Pr[jj]<999)
          {
            sprintf(SS,"%d",Pr[jj]);
/** **/  /*printf("GV2: %s\n",SS);*/
						/*strcat(SS," [label=\"\",color=blue,weight=\n");*/
            strcat(SS," [label=\"\",color=black]\n");
/** **/	/*printf("GV3: %s\n",SS);*/
						strcat(SSSS,SSS);
            strcat(SSSS,SS);
            fprintf(GV,"%s",SSSS);
/** **/ /*printf("GV4:%s\n",SSSS);*/
          }
          if(Pr[jj]==999)
					{
						sprintf(SS,"%d",Pr[jj]);
/** **/ /*printf("GV2: %s\n",SS);*/
						strcat(SS," [label=\"D\",color=black]\n");
/** **/	/*printf("GV3: %s\n",SS);*/
						strcat(SSSS,SSS);
            strcat(SSSS,SS);
            fprintf(GV,"%s",SSSS);
/** **/ /*printf("GV4:%s\n",SSSS);*/
					}
/** **/ /*printf("jj vege: %d\n",jj);*/
        }
      }
    }
  
}
