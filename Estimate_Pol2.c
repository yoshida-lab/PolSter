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
/* SUB ROUTINE OF ESTIMATION */

/* Sparse prediction */
/* Prediction (including missing, exon skipping) */
// pf average

/* State parameter pol2, particle filter generation */
double pol2pfData(int const Par, int n, double pr_mu, double pr_tau, double pr_sigma, char* p){
	printf("check\n");
	
	/* Directory open */
	/* directory name */
	DIR *dp;
	struct dirent *entry;
	char filename[1000];
	//char path[512] = "/home/Kawamura/modeling/NTS_GNPF_Para001/ERR042386_bin400bp_joint_ei_av_chr_sig/";
	//char path[512] = "/mnt/HDPX/ISM/project/ana/modeling/NTS_GNPF_Para003/SRR960177_bin400bp_joint_ei_av_chr_sig_len5cor05_cov01raw/"; // pingu
	//char path[512] = "/work/ana_data/ISM/Modeling/NTS_GNPF_Para003/pf10m/SRR960177_bin400bp_joint_ei_av_chr_sig_len5cor05_cov01raw/test/"; // Luna
	//char path[512] = "/home/yumik/work/modeling/NTS_GNPF_Para003/SRR960177_bin400bp_joint_ei_av_chr_sig_len5cor05_cov01raw_2/"; // ismI
	
	char* path;
	path=p;
	
	puts(p);
	
	int len;
	/* directory open */
	if((dp = opendir(path))==NULL){
		perror("opendir");
		exit(-1);
	}
	len=strlen(path);
	
	// ***** N(3<-100<-10000) CASE GENE DATA ***** //
	
	int g,i,j,k,m,q,r,nm;
	
	// * EI indicator, reduced data EI indicator read
	// exon:1,2,3,.. intron:1,2,3,..
	/* file name */
	sprintf(filename, "%d_EI.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	FILE *fp;
	if((fp=fopen(path,"r"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	// Parameter
	int line;
	line=0;
	char c1[30];
	char c2[30]="E";
	int cmp;
	int INnum; // intron number
	INnum=0;
	int EXnum; // exon number
	EXnum=0;
	int icnt;
	icnt=0;
	int ecnt;
	ecnt=0;
	int Eterm[10000]; // exon number -> exon terminal pos
	int Is[10000]; // intron number -> intron start pos
	int Iterm[10000]; // intron number -> intron terminal pos
	int eiIndic[10000]; // pos -> exon:1, intron:0
	int eNum[10000]; // pos -> exon number, intron case -> 0
	int iNum[10000]; // pos -> intron number, exon case -> 0
	int ex_num; // exon intron number

	while(fgets(c1,30,fp)!=NULL){
		strtok(c1, "\n\0"); 
		cmp = strcmp(c1, c2);
		line++;
		// exon
		if(cmp == 0){
			icnt=0;
			eiIndic[line]=1;
			iNum[line]=0;
			ecnt++;
			if(ecnt==1){
				EXnum++;
			}
			eNum[line]=EXnum;
			Eterm[EXnum]=line;
		}
		// intron
		else{
			ecnt=0;
			eiIndic[line]=0;
			eNum[line]=0;
			icnt++;
			if(icnt==1){
				INnum++;
				Is[INnum]=line;
			}
			iNum[line]=INnum;
			Iterm[INnum]=line;
		}
		//printf("%d\t%d\t%d\t%d\n", line, eiIndic[line], eNum[line], iNum[line]);
		//printf("%d\n", cmp);
	}
	fclose(fp);
	// exon intron number
	ex_num=EXnum;
	
	
	// * Total RNA-seq read data read
	/* file name */
	sprintf(filename, "%d_read.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp=fopen(path,"r"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	// Parameter
	double n4;
	line=0;
	double read[10000];
	int ret;
	while((ret=fscanf(fp,"%lf",&n4))!=EOF){
		line++;
		read[line]=n4;
		//printf("%lf\n",n4);
		//printf("%lf\n", read[line]);
	}
	fclose(fp);
	
	
	/* Write file open */
	/* file name */
	FILE *fp4;
	sprintf(filename, "%d_EI_in.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp4=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	// EI indicator file write
	for(i=1;i<=line;i++){// gene length loop
		fprintf(fp4, "%d\t%d\t%d\t%d\n", i, eiIndic[i], eNum[i], iNum[i]);
	}// (gene len loop)
	fclose( fp4 );
	
	
	
	// ***** STATE PARAMETER ESTIMATION ***** //
	// ***** POL2PF SIMULATION ***** //
	// ***** SNPF SIMULATION ***** //
	
	// --- noise ---> when noise parameter
	double mu,sigma;
	double eps;
	double sigma_eps;
	double a;
	
	mu=sqrt(pr_sigma); ////// optimize, INPUT
	sigma=pr_sigma;
	sigma_eps=pr_tau; ///// optimize, INPUT
	printf("%lf\n",sigma_eps);
	a=0.0;
	// --- noise ---
	
	
	// Parameter
	double**pol2pf=ALLOC_MATRIX(double,row,column);
	double**pol2pfNew=ALLOC_MATRIX(double,row,column); // for replacement
	int**Snpf=ALLOC_MATRIX(int,row,column); // splicing pos
	int**SnpfNew=ALLOC_MATRIX(int,row,column); // splicing pos, fpr replacement
	double**lambdaPF=ALLOC_MATRIX(double,row,column);
	double**lambdaPFtmp=ALLOC_MATRIX(double,row,column);
	double**lambdaPFnew=ALLOC_MATRIX(double,row,column); // for replacement
	int**nums=ALLOC_MATRIX(int,row,column);
	double**w=ALLOC_MATRIX(double,row,column);
	double**wt=ALLOC_MATRIX(double,row,column);
	
	int**flagEI=ALLOC_MATRIX(int,row,column);
	int**flagEItrue=ALLOC_MATRIX(int,row,column);
	int**eiNumEItrue=ALLOC_MATRIX(int,row,column);
	
	double sum_w[10000];
	int PS; // gene length
	PS=line; // for reduced data
	int kk; // exon number or intron number
	int SisP; // Sn intron start pos
	int SitP; // Sn intron term pos
	int SetP; // Sn exon term pos
	//int flagEI; // EI True -> 1, False (cannonical) -> 0
	//int eiNumEItrue; // when EI True, EI number
	//int flagEItrue; //*// When EI true, judge IEI
	
	kk=ex_num; // exon number or intron number
	SisP=PS;
	SitP=PS;
	SetP=PS;
	//flagEI=0;
	//flagEItrue=0; //*//
	//eiNumEItrue=0;
	double prob_sum; // prob normalization
	prob_sum=0.0; // initialization
	double prob_cumu; // prob cumulation
	prob_cumu=0.0; // initialization
	double prob_uni2;
	double prbEql; // equal probability
	
	int ex_ill[1000]; // illegal exon, exon number
	int d; // pos at presence exon number for illegal exon
	int nn; // number of illegal exon
	int flagEill; // flag for illegal exon
	flagEill=0;
	double readEterm;
	double readEsn;
	int flagPrb; // not prob select
	flagPrb=0;
	
	// i=PS generation, initialization
	for(m=0;m<Par;m++){
		//pol2pf[PS][m] = 0.01; // 0.01, initialization
		pol2pf[PS][m] = pr_mu; // 0.01, initialization input
		w[PS][m] = 0.01; // norm dist
		lambdaPF[PS][m] = pol2pf[PS][m];
		lambdaPFtmp[PS][m] = pol2pf[PS][m];
		Snpf[PS][m] = PS;
		flagEI[PS][m] = 0;
		flagEItrue[PS][m] = 0;
		eiNumEItrue[PS][m] = 0;
	}
	printf("%lf\n",pr_mu);
	
	// initialization
	for(i=PS-1;i>0;i--){
		for(m=0;m<Par;m++){
			Snpf[i][m] = PS;
		}
	}
	
	// initialization
	for(i=1;i<=PS;i++){
		sum_w[i] = 0.0;
	}
	
	// * Pol2pf, Snpf, lambda pf generation
	for(i=PS-1;i>0;i--){// i(position) loop
		for(m=0;m<Par;m++){// m(pf) loop
			// -- Xn: pol2 density particle
			eps = rand_normal(0.0,sigma_eps);
			//eps = rand_normal(0.0,sqrt(0.1));
			a = log(pol2pf[i+1][m]) + eps;
			pol2pf[i][m] = exp(a);
						
			// -- Sn: splicing pos particle
			if(eiIndic[i]==1){// exon loop
				kk=eNum[i];
				//Snpf[i][m]=PS; // allmost cannonical
				
				SisP=Is[kk];
				SitP=Iterm[kk];
				SetP=Eterm[(kk+1)];
				
				Snpf[i][m]=Snpf[SetP][m]; // initialization
				
				prob_uni2 = genrand_real2();
				
				flagEI[i][m]=flagEI[i+1][m];
				flagEItrue[i][m]=flagEItrue[i+1][m];
				eiNumEItrue[i][m]=eiNumEItrue[i+1][m];
				
				if(flagEI[i+1][m]==0){// EI False, almost cannonical
					
					// ** equal probability is not including illegal exon skip pattern
					// ** prob sakusei
					// judge illegal exon
					nn=0;
					d=Eterm[kk];
					readEterm=read[d];
					for(g=(kk+1);g<=ex_num;g++){
						j=Eterm[g];
						readEsn=read[j];
						if(readEterm<readEsn){
							nn++;
							ex_ill[nn]=g;
						}
					}
					
					prbEql=0.0;
					prbEql=1.0/(double) (ex_num-kk-nn+1);
					
					// selection
					prob_cumu=0.0;
					for(g=kk;g<=ex_num;g++){
						flagEill=0;
						if(nn>0){
							for(q=1;q<=nn;q++){
								r=ex_ill[q];
								if(r==g){
									flagEill++;
								}
							}
						}
						
						j=Eterm[g];
						if(flagEill==0){
							prob_cumu=prob_cumu+prbEql;
						}
						else{
							prob_cumu=prob_cumu;
						}
						
						flagPrb=0;
						if(prob_cumu>prob_uni2){
							if(g==kk){// intron kk
								Snpf[i][m]=Snpf[i+1][m];
							}
							else{// exon kk+1 ~ ex_num
								Snpf[i][m]=Snpf[j][m];
							}
							flagPrb++;
							break;
						}
					}
					// not prob select
					if(flagPrb==0){
						Snpf[i][m]=Snpf[i+1][m];
					}
					
				}
				else{// EI True, selection from kk to kk_thr (before Enum (Inum) of EI True)
					
					// ** equal probability is not including illegal exon skip pattern
					// ** prob sakusei
					// judge illegal exon
					nn=0;
					d=Eterm[kk];
					readEterm=read[d];
					for(g=(kk+1);g<=ex_num;g++){
						j=Eterm[g];
						readEsn=read[j];
						if(readEterm<readEsn){
							nn++;
							ex_ill[nn]=g;
						}
					}
					prbEql=0.0;
					prbEql=1.0/(double) (eiNumEItrue[i+1][m]-kk-nn+1);
					
					// selection
					prob_cumu=0.0;
					for(g=kk;g<=eiNumEItrue[i+1][m];g++){
						flagEill=0;
						if(nn>0){
							for(q=1;q<=nn;q++){
								r=ex_ill[q];
								if(r==g){
									flagEill++;
								}
							}
						}
						
						j=Eterm[g];
						if(flagEill==0){
							prob_cumu=prob_cumu+prbEql;
						}
						else{
							prob_cumu=prob_cumu;
						}
						
						flagPrb=0;
						if(prob_cumu>prob_uni2){
							if(g==kk){// intron kk
								Snpf[i][m]=Snpf[i+1][m];
								flagEItrue[i][m]=1; //*//
							}
							else{// exon kk+1 ~ ex_num
								Snpf[i][m]=Snpf[j][m];
							}
							flagPrb++;
							break;
						}
					}
					// not prob select
					if(flagPrb==0){
						Snpf[i][m]=Snpf[i+1][m];
					}
				}
			}// (exon loop)
			else{// intron loop
				kk=iNum[i];
				
				SisP=Is[kk];
				SitP=Iterm[kk];
				SetP=Eterm[(kk+1)];
				Snpf[i][m]=SitP; // almost intron terminal, cannonical
				
				prob_uni2 = genrand_real2();
				
				flagEI[i][m]=flagEI[i+1][m];
				flagEItrue[i][m]=flagEItrue[i+1][m];
				eiNumEItrue[i][m]=eiNumEItrue[i+1][m];
				
				if(i==SitP){// Intron terminal
					
					if(flagEItrue[i+1][m]==1){//*// IEI case or not
						if(prob_uni2>0.001){
							Snpf[i][m]=Snpf[(i+1)][m]; //*// IEI case
							flagEI[i][m]=0;
						}
						else{
							Snpf[i][m]=i;
						}
						flagEItrue[i][m]=0; // ?????
					}
					else{
						if(prob_uni2>0.5){// 0.5
							Snpf[i][m]=i;
						}
						else{
							Snpf[i][m]=Snpf[(i+1)][m];
						}
					}
					
					
					// EI True or False judge
					if(Snpf[SetP][m]!=PS){// judge
						if(Snpf[SitP][m]!=Snpf[SetP][m]){//
							flagEI[i][m]=1; // EI True
							eiNumEItrue[i][m]=kk+1; // when EI True, kk threshold, kk+1 ?
						}
					}
					// (judge)
				}
				else{// Intron not terminal
					if(prob_uni2>0.005){// 0.001, * much intron read ??
						Snpf[i][m]=Snpf[(i+1)][m];
					}
					else{
						Snpf[i][m]=i; // maybe RSS
					}
				}
				
			}// (intron loop)
			
			// lambda pf generation
			// lambda exon, intron
			if(eiIndic[i] == 1){// exon loop
				for(j=Snpf[i][m];j>=i;j--){
					if(j==Snpf[i][m]){
						lambdaPFtmp[j][m] = pol2pf[j][m];
					}
					else{
						lambdaPFtmp[j][m] = lambdaPFtmp[j+1][m] + pol2pf[j][m];
					}
				}
				lambdaPF[i][m]=lambdaPFtmp[i][m];
				
			}// (exon loop)
			else{// intron loop
				if(iNum[i] == iNum[i+1]){
					for(j=Snpf[i][m];j>=i;j--){
						if(j==Snpf[i][m]){
							lambdaPFtmp[j][m] = pol2pf[j][m];
						}
						else{
							lambdaPFtmp[j][m] = lambdaPFtmp[j+1][m] + pol2pf[j][m];
						}
					}
					lambdaPF[i][m]=lambdaPFtmp[i][m];
				}
				else{
					lambdaPF[i][m] = pol2pf[i][m];
				}
			}// (intron loop)
			//lambdaPF[i][m] += pol2pf[j][m];
			//printf( "%lf\t%d\t%d\n", lambdaPF[i][m], i, m ); // kakunin
			
			// weight 重み
			if(lambdaPF[i][m] > 0.0){
				wt[i][m] = read[i] / lambdaPF[i][m];
			}
			else{
				lambdaPF[i][m]=0.01;
				wt[i][m] = read[i] / lambdaPF[i][m];
			}
			
			///// likelihood
			if(read[i] > 0.0){
				w[i][m] = (1.0/(sqrt(2.0*M_PI)*sigma*wt[i][m]))*pow(M_E,(-(pow((log(wt[i][m])-mu),2.0))/(2.0*pow(sigma,2.0))));
				// for weight normalized
				sum_w[i] += w[i][m];
			}
			
		}// (m(pf) loop)
		
		// weight normalized 正規化
		if(read[i] > 0.0){
			for(m=0;m<Par;m++){
				w[i][m] = ( w[i][m] / sum_w[i] );
			}
			// SRS resampling algorithm (RSR)
			double u =genrand_res53()/(double)Par;
			
			for(m=0;m<Par;m++){
				nums[i][m]=floor((w[i][m]-u)*Par)+1;
				u+=(nums[i][m]/(double)Par)-w[i][m];
			}
			
			// replacement index, pol2pfNew, SnpfNew, lambdaPF
			for(q=i;q<=PS;q++){
				j=0;
				for(m=0;m<Par;m++){
					nm = nums[i][m];
					for(r=0;r<nm;r++){
						if(j<Par){
							pol2pfNew[q][j] = pol2pf[q][m];
							SnpfNew[q][j] = Snpf[q][m]; //
							lambdaPFnew[q][j] = lambdaPF[q][m];
						 	j++;
						}
					}
				}
			}
			
			// replacement, modoshi, pol2, Snpf
			for(q=i;q<=PS;q++){
				for(m=0;m<Par;m++){
					pol2pf[q][m] = pol2pfNew[q][m];
					Snpf[q][m] = SnpfNew[q][m];
					lambdaPF[q][m] = lambdaPFnew[q][m];
					//printf( "%d\t%d\t%d\n", Snpf[q][m], q, m ); // kakunin
				}
			}
		}
	}// (i(position) loop)
	
	/* Write file open */
	/* file name */
	FILE *fp3;
	sprintf(filename, "%d_pol2pf.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp3=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp14;
	sprintf(filename, "%d_pol2pf_av.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp14=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp15;
	sprintf(filename, "%d_pol2pf_muP2sigma.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp15=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp16;
	sprintf(filename, "%d_pol2pf_muM2sigma.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp16=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	// * Pol2pf value
	// after replacement, pol2pf average, av + 2 sigma, av - 2sigma
	double sum_pf;
	double av_pol2pf[10001];
	double sigma_pol2pf[101];
	double sum_v; // for var cal, pol2, snpf, lambdaPF
	for(i=1;i<=PS;i++){// gene len loop
		// average cal
		sum_pf=0;
		q=0;
		for(m=0;m<Par;m++){// pf loop
			sum_pf += pol2pf[i][m];
			q++;
			fprintf( fp3, "%lf\t%d\t%d\n", pol2pf[i][m],i,m);
		}// (pf loop)
		av_pol2pf[i] = sum_pf / q;
		fprintf( fp14, "%d\t%lf\n", i, av_pol2pf[i]);
		
		// sigma cal
		sum_v=0.0;
		q=0;
		for (m=0;m<Par;m++) {
			sum_v += (pol2pf[i][m] - av_pol2pf[i])*(pol2pf[i][m] - av_pol2pf[i]);
			q++;
		}// (pf loop)
		
		sigma_pol2pf[i] = sqrt(sum_v / q);
		fprintf( fp15, "%d\t%lf\n", i, (av_pol2pf[i] + 2 * sigma_pol2pf[i]) );
		fprintf( fp16, "%d\t%lf\n", i, (av_pol2pf[i] - 2 * sigma_pol2pf[i]) );
		
	}// (gene len loop)
	fclose( fp3 );
	fclose( fp14 );
	fclose( fp15 );
	fclose( fp16 );
	
	
	/* Write file open */
	/* file name */
	FILE *fp13;
	sprintf(filename, "%d_Snpf.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp13=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp12;
	sprintf(filename, "%d_Snpf_av.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp12=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp17;
	sprintf(filename, "%d_Snpf_muP2sigma.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp17=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp18;
	sprintf(filename, "%d_Snpf_muM2sigma.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp18=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	
	// * Snpf value
	// after replacement, Snpf average, av + 2 sigma, av - 2sigma
	int sum_Snpf;
	double av_Snpf[10001];
	double sigma_Snpf[10001];
	for(i=1;i<=PS;i++){// gene len loop
		sum_Snpf=0;
		q=0;
		for(m=0;m<Par;m++){// pf loop
			sum_Snpf += Snpf[i][m];
			q++;
			fprintf( fp13, "%d\t%d\t%d\n", Snpf[i][m],i,m);
		}// (pf loop)
		av_Snpf[i] = (double)sum_Snpf / (double)q;
		fprintf( fp12, "%d\t%lf\n", i, av_Snpf[i]);
		
		// sigma cal
		sum_v=0.0;
		q=0;
		for (m=0;m<Par;m++) {
			sum_v += (Snpf[i][m] - av_Snpf[i])*(Snpf[i][m] - av_Snpf[i]);
			q++;
		}// (pf loop)
		
		sigma_Snpf[i] = sqrt(sum_v / q);
		fprintf( fp17, "%d\t%lf\n", i, (av_Snpf[i] + 2 * sigma_Snpf[i]) );
		fprintf( fp18, "%d\t%lf\n", i, (av_Snpf[i] - 2 * sigma_Snpf[i]) );
		
	}// (gene len loop)
	fclose( fp13 );
	fclose( fp12 );
	fclose( fp17 );
	fclose( fp18 );
	
	
	/* Write file open */
	/* file name */
	FILE *fp5;
	sprintf(filename, "%d_lambdaPF.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp5=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp6;
	sprintf(filename, "%d_lambdaPF_av.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp6=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp19;
	sprintf(filename, "%d_lambdaPF_muP2sigma.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp19=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	/* Write file open */
	/* file name */
	FILE *fp20;
	sprintf(filename, "%d_lambdaPF_muM2sigma.txt", n);
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	if((fp20=fopen(path,"w"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	// * lambda pf value
	// after replacement, lambda pf average, average, av + 2 sigma, av - 2sigma
	double sum_lmd;
	double av_lmdpf[100001];
	double sigma_lmdpf[100001];
	for(i=1;i<=PS;i++){// gene len loop
		sum_pf=0;
		q=0;
		for(m=0;m<Par;m++){// pf loop
			sum_pf += lambdaPF[i][m];
			q++;
			fprintf( fp5, "%lf\t%d\t%d\n", lambdaPF[i][m],i,m);
		}// (pf loop)
		av_lmdpf[i] = sum_pf / q;
		fprintf( fp6, "%d\t%lf\n", i, av_lmdpf[i]);
		
		// sigma cal
		sum_v=0.0;
		q=0;
		for (m=0;m<Par;m++) {
			sum_v += (lambdaPF[i][m] - av_lmdpf[i])*(lambdaPF[i][m] - av_lmdpf[i]);
			q++;
		}// (pf loop)
		
		sigma_lmdpf[i] = sqrt(sum_v / q);
		
		fprintf( fp19, "%d\t%lf\n", i, (av_lmdpf[i] + 2 * sigma_lmdpf[i]) );
		fprintf( fp20, "%d\t%lf\n", i, (av_lmdpf[i] - 2 * sigma_lmdpf[i]) );
		
	}// (gene len loop)
	fclose( fp5 );
	fclose( fp6 );
	fclose( fp19 );
	fclose( fp20 );
	closedir( dp );
	
	*(path+len)='\0'; // ファイル名を付け加える前にpathを最初の状態に戻す
	
	FREE(lambdaPF);
	FREE(lambdaPFnew);
	FREE(lambdaPFtmp);
	FREE(pol2pf);
	FREE(pol2pfNew);
	FREE(Snpf);
	FREE(SnpfNew);
	FREE(nums);
	FREE(wt);
	FREE(w);
	
	return 0;
}
/* MAIN */
int main(int argc, char* argv[]){
	srand((unsigned)time(NULL));
	/* State parameter pol2, Sn estimation, particle filter generation */
	
	/* FILE READ */
	
	/* Directory open */
	/* directory name */
	DIR *dp;
	struct dirent *entry;
	char filename[1000];
	//char path[512] = "/home/Kawamura/modeling/NTS_GNPF_Para001/ERR042386_bin400bp_joint_ei_av_chr_sig/";
	//char path[512] = "/mnt/HDPX/ISM/project/ana/modeling/NTS_GNPF_Para003/SRR960177_bin400bp_joint_ei_av_chr_sig_len5cor05_cov01raw/"; // pingu
	//char path[512] = "/work/ana_data/ISM/Modeling/NTS_GNPF_Para003/pf10m/SRR960177_bin400bp_joint_ei_av_chr_sig_len5cor05_cov01raw/test/"; // Luna
	//char path[512] = "/home/yumik/work/modeling/NTS_GNPF_Para003/SRR960177_bin400bp_joint_ei_av_chr_sig_len5cor05_cov01raw_2/"; // ismI
	
	char* path;
	path=argv[1];
	
	/* directory open */
	int len;
	if((dp = opendir(path))==NULL){
		perror("opendir");
		exit(-1);
	}
	len=strlen(path);
	
	
	// Parameter
	int PS; // gene length
	int line; // gene length count
	int n;
	int g,i,j,k,m,q,r,nm;
	int temp;
	double prob_uni;
	double prob_uni2;
	double prob_sum; // prob normalization
	prob_sum=0.0;
	double prob_cumu; // prob cumulation
	prob_cumu=0.0;
	double prbNrm[10000]; // tashite 1 ninaru probs
	double prbEql; // equal probability
	int eiIndic[10000]; // pos -> exon:1, intron:0
	int eNum[10000]; // pos -> exon number, intron case -> 0
	int iNum[10000]; // pos -> intron number, exon case -> 0
	int ex_num; // exon intron number
	
	
	// * Hyper parameter mu, tau, sigma INPUT read
	/* File open */
	/* file name */
	sprintf(filename, "SRR960177_pickup_gene_num_joint_ei_av_mu_tau_sigma.txt");
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	FILE *fp;
	if((fp=fopen(path,"r"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	
	// Data
	int n10, n11;
	double n12, n13, n14;
	double para_mu[10000]; // pol2pf inititalization
	double para_tau[10000]; // pol2pf noise eps
	double para_sigma[10000]; // weight cal noise sigma
	int ret;
	while((ret=fscanf(fp,"%d\t%d\t%lf\t%lf\t%lf",&n10,&n11,&n12,&n13,&n14))!=EOF){
		para_mu[n11]=n12;
		para_tau[n11]=n13;
		para_sigma[n11]=n14;
		printf("%d\t%d\t%lf\t%lf\t%lf\n", n10,n11,n12,n13,n14);
	}
	fclose(fp);
	
	
	// * Partcle number for PF INPUT read
	/* File open */
	/* file name */
	sprintf(filename, "Particle_num.txt");
	*(path+len)='\0';
	strcat(path, filename);
	/* file open */
	FILE *fp2;
	if((fp2=fopen(path,"r"))==NULL){
		printf("%s can't open.\n",filename);
	}
	else {
		printf("%s open.\n",filename);
	}
	*(path+len)='\0'; // ファイル名を付け加える前にpathを最初の状態に戻す
	
	// Data
	int Par; // particle number for PF
	//Par = 100000; // 100000
	int ret2;
	while((ret2=fscanf(fp2,"%d",&Par))!=EOF){
		printf("Particle num: %d\n", Par);
	}
	fclose(fp2);
	closedir( dp );
	
	// ***** Data start ***** //
	double pr_mu; // pol2pf inititalization
	double pr_tau; // pol2pf noise eps
	double pr_sigma; // weight cal noise sigma
	
	// * n => num of genes
	for(n=10;n<=19;n++){// n loop
		printf("check\n");
		pr_mu = para_mu[n];
		pr_tau = para_tau[n];
		pr_sigma = para_sigma[n];
		if(pr_mu==0.0){
			printf("check\n");
		}
		else{
			pol2pfData(Par,n,pr_mu,pr_tau,pr_sigma,path); // -> sub routine of estimation
		}
	}// (n loop)
	
	return 0;
}
