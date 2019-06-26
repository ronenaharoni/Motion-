#include <math.h>
#include <stdlib.h>
#include "MotionHeader.h"
#include "polyfit.h"

int MotionTracking(float* Mscan[],SysParams_Struct* SysParams){
	int i,k,t,m;
	int DistThr= 5;
	float Centoird1,Centoird2, Noise_dB;
	float SumOfSecondDiff,MeanOfDiff,FirstDiff[SysParams->Nscans-1],SecondDiff[SysParams->Nscans-2];
	float SumForVar,NoiseEstimation,PowerOfMscan,maxEnergyPerTime;
	//	float energyPerTime_dB[SysParams->Nbins][SysParams->SpectrogramTimeBinsSingleMotion];
	double 	ValuesForInterp[SysParams->SpectrogramTimeBinsSingleMotion];
	double PointsForInterp[SysParams->SpectrogramTimeBinsSingleMotion];
	int NumPointsForInterp=0;
	double p[2];//1st order polynom coefficients
	float   NextRange,RstartWanted;
	float* energyPerTime_dB[SysParams->Nbins];//,*Mscan_old[SysParams.Nscans];
	for (i=0; i<SysParams->Nbins; i++){
		energyPerTime_dB[i] = (float *)malloc(SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	}

	//	int FirstBinAboveThresholdPerTime[SysParams->SpectrogramTimeBinsSingleMotion];
	//	memset(FirstBinAboveThresholdPerTime,0,(SysParams->SpectrogramTimeBinsSingleMotion)*sizeof(int));

	for(i=0;i<SysParams->Nbins;i++){//calculate NoiseEstimation=var(diff(Signal_k,2)) for each bin
		SumOfSecondDiff=0;
		for(k=0;k<SysParams->Nscans-1;k++){//1st Diff
			FirstDiff[k]=(Mscan[k+1][i]-Mscan[k][i]);
		}
		for(k=0;k<SysParams->Nscans-2;k++){//2nd Diff
			SecondDiff[k]=(FirstDiff[k+1]-FirstDiff[k]);
			SumOfSecondDiff+=SecondDiff[k];
			//			printf("%d %lf\n",k,SecondDiff[k]);
		}

		MeanOfDiff=SumOfSecondDiff/(SysParams->Nscans-2);//Calculate mean
		SumForVar=0;
		for(k=0;k<SysParams->Nscans-2;k++){//VAR=SUM((x-mu)^2)/N-1,  N-1 because it's unbiased estimator
			SumForVar+=pow((SecondDiff[k]-MeanOfDiff),2);
		}
		NoiseEstimation=SumForVar/((SysParams->Nscans-2)-1);//it's the Variance

		for(t=0;t<SysParams->SpectrogramTimeBinsSingleMotion;t++){//take 30 points of slow and calculate the max energy of them
			maxEnergyPerTime=0;
			for(m=0;m<SysParams->winLength;m++){
				PowerOfMscan=pow(Mscan[m+t*(SysParams->TimeShift)][i],2);
				//				printf("%d %lf\n",m,Mscan[m+t*(SysParams->TimeShift)][i]);

				if(PowerOfMscan > maxEnergyPerTime){
					maxEnergyPerTime=PowerOfMscan;
				}
			}

			energyPerTime_dB[i][t]=10*log10(maxEnergyPerTime/NoiseEstimation);
			//			printf("%d %d    %lf\n",i,t,energyPerTime_dB[i][t]);

		}

	}

	K_means(energyPerTime_dB,SysParams,&Centoird1,&Centoird2);//Calulcate the centroids of the samples

	if(fabs(Centoird1-Centoird2)>DistThr){
		if(Centoird1>Centoird2)
			Noise_dB =Centoird1;
		else
			Noise_dB =Centoird2;
	}
	else
		Noise_dB=INFINITY;//why??



	//1st order interpolation for the curve of the walking

	for(t=0;t<SysParams->SpectrogramTimeBinsSingleMotion;t++){
		for(i=0;i<SysParams->Nbins;i++){
			if(energyPerTime_dB[i][t]>Noise_dB){//find the first Bin which crosses the energy
				ValuesForInterp[NumPointsForInterp]=SysParams->Rbin_m[i];//this is the first rbin_m who cross the threshold put in rbin_m
				PointsForInterp[NumPointsForInterp]=t+1;//save the valid time point for interpolation, +1 to be like matlab
				NumPointsForInterp+=1;
				break;
			}
		}

	}

	//Add this
	if((SysParams->SpectrogramTimeBinsSingleMotion - NumPointsForInterp) > 0.87*SysParams->SpectrogramTimeBinsSingleMotion){//if we don't have enough points -> don't change
		SysParams->Rstart_corrected=SysParams->Rbin_m[0];
		SysParams->Rstop_corrected=SysParams->Rbin_m[SysParams->SpectrogramTimeBinsSingleMotion-1]-0.1;
		return 0;
	}

	polyfit(PointsForInterp,ValuesForInterp,NumPointsForInterp,1,p);//1 si the order MUST TO BE float

	NextRange=p[1]*SysParams->SpectrogramTimeBinsSingleMotion+p[0];//insert the last point in order to detemine the center of the next ring

	RstartWanted=SysParams->Rmin;
	if ((NextRange-0.5*2.5) > RstartWanted)//max(nextRange - 1.25, Rmin)
		RstartWanted=NextRange-0.5*2.5;

	if((RstartWanted + 2.5) > SysParams->Rmax){//if the end of the ring is exit from the Rmax, cut it by the Rmax
		SysParams->Rstart_corrected = SysParams->Rmax - 2.5;
		SysParams->Rstop_corrected=SysParams->Rmax;
	}
	else{
		SysParams->Rstart_corrected = RstartWanted;
		SysParams->Rstop_corrected = RstartWanted + 2.5;
	}


	return 0;
}


int K_means(float* energyPerTime_dB[],SysParams_Struct* SysParams,float* ChosenCentoird1,float* ChosenCentoird2){
	float SumForCenteroid1, SumForCenteroid2,Centeroid1,Centeroid2,PrevCenteroid1,PrevCenteroid2;
	 float epsilon=0;
	int i,t,AlgorithmDone=0;
	int NumOfSamplesCentroid1,NumOfSamplesCentroid2,NumOfIter=0;
	SumForCenteroid1=0;
	SumForCenteroid2=0;
	int Half_K_mean_TotalSamples=SysParams->SpectrogramTimeBinsSingleMotion * SysParams->Nbins/2;
	//find initial centroids by deviding to two groups
	for(i=0;i<SysParams->Nbins;i++){
		for(t=0;t<SysParams->SpectrogramTimeBinsSingleMotion;t++){
			if(t<SysParams->SpectrogramTimeBinsSingleMotion/2)//first half
				SumForCenteroid1+=energyPerTime_dB[i][t];
			else//second half
				SumForCenteroid2+=energyPerTime_dB[i][t];
		}
	}
	Centeroid1=SumForCenteroid1/(Half_K_mean_TotalSamples);
	Centeroid2=SumForCenteroid2/(Half_K_mean_TotalSamples);


	while(AlgorithmDone==0){
		SumForCenteroid1=0;
		SumForCenteroid2=0;
		NumOfSamplesCentroid1=0;
		NumOfSamplesCentroid2=0;

		for(i=0;i<SysParams->Nbins;i++){
			for(t=0;t<SysParams->SpectrogramTimeBinsSingleMotion;t++){

//				printf("%d %d    %lf\n",i,t,energyPerTime_dB[i][t]);

				if(fabs(energyPerTime_dB[i][t] - Centeroid1) < fabs(energyPerTime_dB[i][t] - Centeroid2)){//if the sample is closer to centroid 1
					SumForCenteroid1+=energyPerTime_dB[i][t];
					NumOfSamplesCentroid1+=1;
				}
				else{//the sample is closer to centroid2
					SumForCenteroid2+=energyPerTime_dB[i][t];
					NumOfSamplesCentroid2+=1;
				}
			}
		}
			//new centroids:
			PrevCenteroid1=Centeroid1;
			PrevCenteroid2=Centeroid2;

			Centeroid1=SumForCenteroid1/NumOfSamplesCentroid1;
			Centeroid2=SumForCenteroid2/NumOfSamplesCentroid2;


		//Stop condition
		if ((fabs(Centeroid1-PrevCenteroid1) <= epsilon) && (fabs(Centeroid2-PrevCenteroid2) <= epsilon)){
			*ChosenCentoird1=Centeroid1;
			*ChosenCentoird2=Centeroid2;
			AlgorithmDone=1;
		}
		NumOfIter++;
	}
	return 0;

}
