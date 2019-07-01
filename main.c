#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "MotionHeader.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "test.h"
#include "seatest.h"


//struct timeval tpStart , tpStop ;
//    struct timeval tp ;
//	float f1 = 0 ;
//	gettimeofday(&tp,0);

int main( void ){
	//	FILE *MscanFile = fopen("Mscan_withoutMF202.csv", "r");
	FILE *MscanFile1 = fopen("/home/debian/Records/Mscan_withoutMF700_1.csv", "r");
	FILE *MscanFile2 = fopen("/home/debian/Records/Mscan_withoutMF700_2.csv", "r");
	FILE *MscanFile3 = fopen("/home/debian/Records/Mscan_withoutMF700_3.csv", "r");

	SysParams_Struct SysParams;

	SysParams.Rmin=1;//need to change this from the configuration XML!!!!
	SysParams.Ring=2.62572418212891; //the ring width in meters; FOR DEMO ONLY maybe 2.5??? why 2.6
	SysParams.Rmax=3.5; //need to change this from the configuration XML!!!!

	SysParams.FirstTimeMotion=1;//1 if it's the first 2 frames
	SysParams.Nbins=288;
	SysParams.Fs=250;
	SysParams.Nscans=400; // 400 for 1.6 sec record, 800 in the single 3.2 sec record
	SysParams.NumOfClasses=4;
	SysParams.RangeWide=120;
	SysParams.winLength=30;
	SysParams.N_Overlap=25;
	SysParams.TimeShift=SysParams.winLength-SysParams.N_Overlap;//=5
	SysParams.Fbins=100;
	SysParams.noiseFreqBins=4;
	SysParams.Fmin_bin=2;
	SysParams.MinFreq=2.5;//in Hz
	SysParams.minEventDuration=30;
	SysParams.DFTLengthForPSD=floor(SysParams.Nscans/2)+1;
	SysParams.noiseThresh=5;
	SysParams.DFTLengthForSpectrogram=2*SysParams.Fbins-1;//=199
	SysParams.GapLength=14;
	SysParams.SpectrogramFreqBins=SysParams.DFTLengthForSpectrogram/2+1;//=100
	SysParams.SpectrogramTimeBinsSingleMotion = floor((SysParams.Nscans-SysParams.N_Overlap)/(SysParams.winLength - SysParams.N_Overlap)); // total time bins of the spectrogram=155
	SysParams.SpectrogramTimeBinsTwoMotions=(2*SysParams.SpectrogramTimeBinsSingleMotion+SysParams.GapLength);//=2*75+14=164
	SysParams.SpectrogramTimeBinsThreeMotions=(3*SysParams.SpectrogramTimeBinsSingleMotion+2*SysParams.GapLength);//=3*75+2*14=253
	SysParams.SpectrogramFreqBinsHilbert=2*SysParams.Fbins;//=200
	SysParams.DFTLengthForSpectrogramHilbert=2*SysParams.Fbins;//=200
	SysParams.MedianValue=20;//for MedianFilter
	SysParams.AvgValue=5;//for AvgFilter
	SysParams.topMaxFreq=5;
	SysParams.truncate=0;// Computes medians of smaller segments as it reaches the signal edges in Pxx2 (regular) Spectrogram
	SysParams.truncateHilbert=1;// Computes medians of smaller segments as it reaches the signal edges in Hilbert Spectrogram

	//Hamming window coefficients
	SysParams.Hamming[0]=0.08; SysParams.Hamming[1]=0.0907545443733602; SysParams.Hamming[2]=0.122515306951360;
	SysParams.Hamming[3]=0.173797189775404; SysParams.Hamming[4]=0.242202309000359; SysParams.Hamming[5]=0.324532117278097;
	SysParams.Hamming[6]=0.416936964276558;SysParams.Hamming[7]=0.515096102050708; SysParams.Hamming[8]=0.614419718414272;
	SysParams.Hamming[9]=0.710263551456361; SysParams.Hamming[10]=0.798146050066696;SysParams.Hamming[11]=0.873957926284640;
	SysParams.Hamming[12]=0.934154301037091; SysParams.Hamming[13]=0.975920458744089; SysParams.Hamming[14]=0.997303460291006;
	SysParams.Hamming[15]=SysParams.Hamming[14]; SysParams.Hamming[16]=SysParams.Hamming[13]; SysParams.Hamming[17]=SysParams.Hamming[12]; SysParams.Hamming[18]=SysParams.Hamming[11];
	SysParams.Hamming[19]=SysParams.Hamming[10]; SysParams.Hamming[20]=SysParams.Hamming[9]; SysParams.Hamming[21]=SysParams.Hamming[8]; SysParams.Hamming[22]=SysParams.Hamming[7];
	SysParams.Hamming[23]=SysParams.Hamming[6]; SysParams.Hamming[24]=SysParams.Hamming[5]; SysParams.Hamming[25]=SysParams.Hamming[4]; SysParams.Hamming[26]=SysParams.Hamming[3];
	SysParams.Hamming[27]=SysParams.Hamming[2]; SysParams.Hamming[28]=SysParams.Hamming[1]; SysParams.Hamming[29]=SysParams.Hamming[0];
	SysParams.A12_Inverse[0][0]=5.925925925925995e-04;SysParams.A12_Inverse[0][1]=-5.925925925925996e-04;SysParams.A12_Inverse[0][2]=0.004444444444445;SysParams.A12_Inverse[0][3]=0.004444444444444;
	SysParams.A12_Inverse[1][0]=-0.146666666666668;SysParams.A12_Inverse[1][1]=0.146666666666668;SysParams.A12_Inverse[1][2]=-1.133333333333347;SysParams.A12_Inverse[1][3]=-1.066666666666678;
	SysParams.A12_Inverse[2][0]=12.0000000000001;SysParams.A12_Inverse[2][1]=-12.0000000000001;SysParams.A12_Inverse[2][2]=96.0000000000011;SysParams.A12_Inverse[2][3]=85.0000000000010;
	SysParams.A12_Inverse[3][0]=-324.000000000004;SysParams.A12_Inverse[3][1]=325.000000000004;SysParams.A12_Inverse[3][2]=-2.700000000000030e+03;SysParams.A12_Inverse[3][3]=-2.250000000000026e+03;
	SysParams.NumSamplesForDerivativeEstimation=5;

	Motion_Struct2 MotionStruct0;


	Edge2_Struct Edge2_0;
	Edge2_Struct Edge2_Plus_0;
	Edge2_Struct Edge2_Minus_0;

	int i,j;
	int NumOfTrees=8,NumOfFeatures=7;
	float delta_R;
	float Mscan_flat[SysParams.Nscans*SysParams.Nbins];
	float* Mscan0[SysParams.Nscans],* Mscan1[SysParams.Nscans],*Mscan2[SysParams.Nscans],* Mscan3[SysParams.Nscans];//,*Mscan_old[SysParams.Nscans];
	for (i=0; i<SysParams.Nscans; i++){
		Mscan0[i] = (float *)malloc(SysParams.Nbins * sizeof(float));
		Mscan1[i] = (float *)malloc(SysParams.Nbins * sizeof(float));
		Mscan2[i] = (float *)malloc(SysParams.Nbins * sizeof(float));
		Mscan3[i] = (float *)malloc(SysParams.Nbins * sizeof(float));

	}
	Tree_Struct Tree[8];
	Tree_Struct* All_Trees[8]={&Tree[0],&Tree[1],&Tree[2],&Tree[3],&Tree[4],&Tree[5],&Tree[6],&Tree[7]};
	RF_Struct RF_Model;
	SVM_Struct SVM_Model;

	/////MotioStruct0 Preparation/////
	SysParams.SpectrogramTimeBins=SysParams.SpectrogramTimeBinsSingleMotion;

	Edge2_0.FiftyPrecent_Filtered=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_0.Fmax=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_0.PrevLastFiftyPrecent=(float *)malloc((SysParams.AvgValue-1) * sizeof(float));

	if(SysParams.truncate==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
		Edge2_0.Peak=(float *)malloc((SysParams.SpectrogramTimeBins) * sizeof(float));
		Edge2_0.Peak_Filtered=(float *)malloc((SysParams.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_0.Peak,0,(SysParams.SpectrogramTimeBins+SysParams.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_0.Peak_Filtered=(float *)malloc((SysParams.SpectrogramTimeBins+SysParams.MedianValue-1) * sizeof(float));
	}

	Edge2_Plus_0.FiftyPrecent_Filtered=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.Fmax=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.PrevLastFiftyPrecent=(float *)malloc((SysParams.AvgValue-1) * sizeof(float));
	Edge2_Plus_0.SumEnergy_Post=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.T1_t=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));

	if(SysParams.truncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_0.Peak_Filtered=(float *)malloc((SysParams.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Plus_0.Peak,0,(SysParams.SpectrogramTimeBins+SysParams.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_0.Peak_Filtered=(float *)malloc((SysParams.SpectrogramTimeBins+SysParams.MedianValue-1) * sizeof(float));
	}

	Edge2_Minus_0.FiftyPrecent_Filtered=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.Fmax=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.PrevLastFiftyPrecent=(float *)malloc((SysParams.AvgValue-1) * sizeof(float));
	Edge2_Minus_0.SumEnergy_Post=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.T1_t=(float *)malloc(SysParams.SpectrogramTimeBins * sizeof(float));

	if(SysParams.truncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_0.Peak_Filtered=(float *)malloc((SysParams.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Minus_0.Peak,0,(SysParams.SpectrogramTimeBins+SysParams.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the initial conditions
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_0.Peak_Filtered=(float *)malloc((SysParams.SpectrogramTimeBins+SysParams.MedianValue-1) * sizeof(float));
	}

	//	if(SysParams.FirstTimeMotion==1){//if it's first time there is no history and therefore set all 0 for the initial state of the AvgFilter
	//		memset(Edge2_0.PrevLastFiftyPrecent,0,(SysParams.AvgValue-1)*sizeof(float));
	//		memset(Edge2_Plus_0.PrevLastFiftyPrecent,0,(SysParams.AvgValue-1)*sizeof(float));
	//		memset(Edge2_Minus_0.PrevLastFiftyPrecent,0,(SysParams.AvgValue-1)*sizeof(float));
	//	}


	MotionStruct0.Edge2=&Edge2_0;
	MotionStruct0.Edge2_Plus=&Edge2_Plus_0;
	MotionStruct0.Edge2_Minus=&Edge2_Minus_0;

	RF_Params_Import(NumOfTrees, NumOfFeatures,All_Trees,&RF_Model);
	SVM_Params_Import(&SVM_Model);

	read_data_from_file(MscanFile1,SysParams.Nscans,SysParams.Nbins,Mscan_flat);//get the MSCAN in flat form
	for ( i=0;i<SysParams.Nscans;i++)
	{//get the  MSCAN from flat
		for ( j=0;j<SysParams.Nbins;j++)
		{
			Mscan0[i][j]=Mscan_flat[i*SysParams.Nbins+j];
			//			Mscan_old[i][j]=Mscan_flat[i*Nbins+j];
		}
	}

	read_data_from_file(MscanFile1,SysParams.Nscans,SysParams.Nbins,Mscan_flat);//get the MSCAN in flat form
	for ( i=0;i<SysParams.Nscans;i++)
	{//get the  MSCAN from flat
		for ( j=0;j<SysParams.Nbins;j++)
		{
			Mscan1[i][j]=Mscan_flat[i*SysParams.Nbins+j];
			//			Mscan_old[i][j]=Mscan_flat[i*Nbins+j];
		}
	}

	//create Rbin_m
	SysParams.Rstart_corrected=2.231486865000000;//for demo only!!! it will be the current Rstart
	SysParams.Rbin_m=(float *)malloc(SysParams.Nbins * sizeof(float));
	delta_R=(SysParams.Ring)/((float)SysParams.Nbins-1);//the delta of the bins inside the ring(=Rstop-Rstart)=2.5m

	for(i=0;i<SysParams.Nbins;i++){
		SysParams.Rbin_m[i]=SysParams.Rstart_corrected+delta_R*i;
//		printf("%d %lf\n",i,SysParams.Rbin_m[i]);
	}



//	gettimeofday(&tpStart,0);
			// end point

	MotionHandler(Mscan0,Mscan1,&SysParams,All_Trees,&SVM_Model,&RF_Model,&MotionStruct0);
//	gettimeofday(&tpStop,0);
//
//	f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000) ;
//			printf (" %f sec\n", f1 );


//	free(Mscan0);NNED TO FREE
//	free(Mscan1);

	////////////////////


	SysParams.FirstTimeMotion=0;//now will be 3 motions


	read_data_from_file(MscanFile2,SysParams.Nscans,SysParams.Nbins,Mscan_flat);//get the MSCAN in flat form
	for ( i=0;i<SysParams.Nscans;i++)
	{//get the  MSCAN from flat
		for ( j=0;j<SysParams.Nbins;j++)
		{
			Mscan2[i][j]=Mscan_flat[i*SysParams.Nbins+j];
		}
	}


	read_data_from_file(MscanFile3,SysParams.Nscans,SysParams.Nbins,Mscan_flat);//get the MSCAN in flat form
	for ( i=0;i<SysParams.Nscans;i++)
	{//get the  MSCAN from flat
		for ( j=0;j<SysParams.Nbins;j++)
		{
			Mscan3[i][j]=Mscan_flat[i*SysParams.Nbins+j];
		}
	}

	MotionHandler(Mscan2,Mscan3,&SysParams,All_Trees,&SVM_Model,&RF_Model,&MotionStruct0);



//	free(Mscan2);
//	free(Mscan3);




	return 0;
}


int read_data_from_file(FILE *fp_read ,int Nscans,int Nbins,float* RespMscan_flat)
{
	int i;


	for (i=0;i<Nscans*Nbins;i++)
	{
		fscanf(fp_read,"%f",RespMscan_flat+i);
		//				printf("%d %lf\n",i,*(RespMscan_flat+i));
	}
	return 0;
}


int RF_Params_Import(int NumOfTrees,int NumOfFeatures, Tree_Struct** All_Trees,RF_Struct* RF_Model){
	FILE *x_mean_File;
	FILE *x_std_File;
	double Value0;
	int i;
	x_mean_File=fopen("/home/debian/RandomForestClassifier/RF_x_mean.csv", "r");
	x_std_File=fopen("/home/debian/RandomForestClassifier/RF_x_std.csv", "r");

	for(i=0;i<NumOfFeatures;i++){//import the mean and std for the normalization
		fscanf(x_mean_File,"%lf",&Value0);
		RF_Model->x_mean[i]=Value0;
		fscanf(x_std_File,"%lf",&Value0);
		RF_Model->x_std[i]=Value0;
	}

	TreesCreator(NumOfTrees, All_Trees);//create 8 trees of the classifier
	return 0;
}

int TreesCreator(int NumOfTrees, Tree_Struct** All_Trees){

	//notice: some variabels are hard coded
	char  Filename_Children[50] ;
	char Filename_CutPoint[50];
	char Filename_CutPredictor[54];
	char Filename_ClassProb[50];
	FILE *Children_File;
	FILE *CutPoint_File;
	FILE *CutPredictor_File;
	FILE *ClassProb_File;
	FILE *Lengths_File;
	int TreeNum,i;
	int Current_Length;//the length of each tree
	float Value0,Value1,Value2,Value3;
	int ValueCutPredictor;
	int Value0Children,Value1Children;
	Lengths_File=fopen("/home/debian/RandomForestClassifier/Lengths.csv", "r");
	All_Trees[0]->TotalTrees=8;
	for(TreeNum=0;TreeNum<All_Trees[0]->TotalTrees;TreeNum++){
		fscanf(Lengths_File,"%d",&Current_Length);
		memset(Filename_Children , 0 , sizeof(Filename_Children));
		memset(Filename_CutPoint , 0 , sizeof(Filename_CutPoint));
		memset(Filename_CutPredictor , 0 , sizeof(Filename_CutPredictor));
		memset(Filename_ClassProb , 0 , sizeof(Filename_ClassProb));

		sprintf(Filename_CutPredictor, "/home/debian/RandomForestClassifier/CutPredictor%d.csv" , TreeNum);
		CutPredictor_File=fopen(Filename_CutPredictor, "r");


		sprintf(Filename_CutPoint, "/home/debian/RandomForestClassifier/CutPoint%d.csv" , TreeNum);
		CutPoint_File=fopen(Filename_CutPoint, "r");


		sprintf(Filename_Children, "/home/debian/RandomForestClassifier/Children%d.csv" , TreeNum);
		Children_File=fopen(Filename_Children, "r");

		sprintf(Filename_ClassProb, "/home/debian/RandomForestClassifier/ClassProb%d.csv" , TreeNum);
		ClassProb_File=fopen(Filename_ClassProb, "r");


		All_Trees[TreeNum]->CutPoint =(float *)malloc(Current_Length * sizeof(float));
		All_Trees[TreeNum]->CutPredictor =(int *)malloc(Current_Length * sizeof(int));
		All_Trees[TreeNum]->Children=(int *)malloc(2 *Current_Length* sizeof(int));
		All_Trees[TreeNum]->ClassProb=(float *)malloc(4 *Current_Length* sizeof(float));

		for(i=0;i<Current_Length;i++){
			fscanf(CutPoint_File,"%f",&Value0);//Value=-1 is for NaN
			//     		printf("%f cut\n",Value0);
			All_Trees[TreeNum]->CutPoint[i]=Value0;
			fscanf(CutPredictor_File,"%d",&ValueCutPredictor);
			All_Trees[TreeNum]->CutPredictor[i]=ValueCutPredictor;//Value=-1 is for NaN
			//			     		printf("%d cutp\n",ValueCutPredictor);
			fscanf(Children_File,"%d,%d",&Value0Children,&Value1Children);

			All_Trees[TreeNum]->Children[2*i]=Value0Children;//children is 2 columns left & right
			All_Trees[TreeNum]->Children[2*i+1]=Value1Children;

			fscanf(ClassProb_File,"%f,%f,%f,%f",&Value0,&Value1,&Value2,&Value3);

			All_Trees[TreeNum]->ClassProb[4*i]=Value0;//classprob is 4 columns because 4 classes
			All_Trees[TreeNum]->ClassProb[4*i+1]=Value1;
			All_Trees[TreeNum]->ClassProb[4*i+2]=Value2;
			All_Trees[TreeNum]->ClassProb[4*i+3]=Value3;
			//			printf("%f,%f,%f,%f class\n",Value0,Value1,Value2,Value3);


		}
	}
	return 0;
}


int SVM_Params_Import(SVM_Struct* SVM_Model){
	FILE *Bias_File;
	FILE *Beta_File;
	FILE *x_mean_File;
	FILE *x_std_File;
	float Value0;
	int i,NumOfFeatures=7;

	Bias_File=fopen("/home/debian/SVMClassifier/SVM_Bias.csv", "r");
	Beta_File=fopen("/home/debian/SVMClassifier/SVM_Beta.csv", "r");
	x_mean_File=fopen("/home/debian/SVMClassifier/SVM_x_mean.csv", "r");
	x_std_File=fopen("/home/debian/SVMClassifier/SVM_x_std.csv", "r");

	for(i=0;i<NumOfFeatures;i++){
		fscanf(Beta_File,"%f",&Value0);
		SVM_Model->Beta[i]=Value0;
		fscanf(x_mean_File,"%f",&Value0);
		SVM_Model->x_mean[i]=Value0;
		fscanf(x_std_File,"%f",&Value0);
		SVM_Model->x_std[i]=Value0;
	}
	fscanf(Bias_File,"%f",&Value0);
	SVM_Model->Bias=Value0;

	return 0;
}
