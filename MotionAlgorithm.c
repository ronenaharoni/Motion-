#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
#define I _Complex_I (0.0f+1.0if)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "MotionHeader.h"
#include <math.h>
#include <gsl/gsl_movstat.h>
#include "polyfit.h"
#include <fftw3.h>

///ADD MSAN 1 AND MSCAN 2 AND AT THE END SAVE MOTION2 AT MOTION0
int MotionHandler(float* Mscan1[], float* Mscan2[], SysParams_Struct* SysParams,
		Tree_Struct** All_Trees, SVM_Struct* SVM_Model, RF_Struct* RF_Model,
		Motion_Struct2* MotionStruct0) {
	int i;
	int y_hat_M;
	float* Mscan_abs_FFT[SysParams->DFTLengthForPSD];
	float* Pxx2_dB[SysParams->SpectrogramFreqBins],
	*Pxx2[SysParams->SpectrogramFreqBins];
	float * Pxx2_Hilbert[SysParams->SpectrogramFreqBinsHilbert];
	float MotionDistribution[SysParams->NumOfClasses]; //=4

	Motion_Struct2 MotionStruct1;
	Motion_Struct2 MotionStruct2;
	Motion_Struct2 UnitedMotionStruct;

	Edge2_Struct Edge2_1;
	Edge2_Struct Edge2_Plus_1;
	Edge2_Struct Edge2_Minus_1;

	Edge2_Struct Edge2_2;
	Edge2_Struct Edge2_Plus_2;
	Edge2_Struct Edge2_Minus_2;

	Edge2_Struct Edge2_United;
	Edge2_Struct Edge2_Plus_United;
	Edge2_Struct Edge2_Minus_United;

	Features_Struct Plus, Minus, Both;
	AllFeatures_Struct FeatureSet;

	Event_Struct EventStruct_United;

	FeatureSet.Plus = &Plus;
	FeatureSet.Minus = &Minus;
	FeatureSet.Both = &Both;

	//////////////////////////Calculation of MotionStruct1//////////////////////////

	//	MotionStruct1.EventStruct=&EventStruct1;
	SysParams->SpectrogramTimeBins = SysParams->SpectrogramTimeBinsSingleMotion; //first there is extraction on single motion (=75)

	for (i = 0; i < SysParams->DFTLengthForPSD; i++) {
		Mscan_abs_FFT[i] = (float *) malloc(SysParams->Nbins * sizeof(float));
	}
	for (i = 0; i < SysParams->SpectrogramFreqBins; i++) {
		Pxx2_dB[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
		Pxx2[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
	}
	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert; i++) {

		Pxx2_Hilbert[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
	}

	//PreProcessing
	MacthedFilter2(Mscan1, SysParams);

	//	SlowProcessing2(Mscan1, SysParams); //remove DC
	//	NotchFilter2(Mscan1, SysParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp

	//Memory Allocation for the curves
	MemoryAllocation(&Edge2_1, &Edge2_2, &Edge2_Plus_1, &Edge2_Minus_1,
			&Edge2_Plus_2, &Edge2_Minus_2, SysParams);

	if (SysParams->FirstTimeMotion == 1) { //if it's first time there is no history and therefore set all 0 for the initial state of the AvgFilter
		memset(Edge2_1.PrevLastFiftyPrecent, 0,
				(SysParams->AvgValue - 1) * sizeof(float));
		memset(Edge2_Plus_1.PrevLastFiftyPrecent, 0,
				(SysParams->AvgValue - 1) * sizeof(float));
		memset(Edge2_Minus_1.PrevLastFiftyPrecent, 0,
				(SysParams->AvgValue - 1) * sizeof(float));
	} else {
		memcpy(Edge2_Plus_1.PrevLastFiftyPrecent,
				MotionStruct0->Edge2_Plus->PrevLastFiftyPrecent,
				(SysParams->AvgValue - 1) * sizeof(float));
		memcpy(Edge2_Minus_1.PrevLastFiftyPrecent,
				MotionStruct0->Edge2_Minus->PrevLastFiftyPrecent,
				(SysParams->AvgValue - 1) * sizeof(float));
		//notice that PrevLastFiftyPrecent is actually the last 50precent of MotionStruct0;
	}

	//Extraction of the curves
	MotionCurveExtraction2(&MotionStruct1, Mscan1, Mscan_abs_FFT, Pxx2_Hilbert,
			Pxx2, Pxx2_dB, &Edge2_1, &Edge2_Plus_1, &Edge2_Minus_1, SysParams);

	//	////////////Send Mscan1 to tracking/////////////////////

	MotionTracking(Mscan1, SysParams);

	///////////FREE MEMORY///////////////

	for (i = 0; i < SysParams->DFTLengthForPSD; i++) {
		free(Mscan_abs_FFT[i]);
	}
	for (i = 0; i < SysParams->SpectrogramFreqBins; i++) {
		free(Pxx2_dB[i]);
		free(Pxx2[i]);
	}
	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert; i++) {
		free(Pxx2_Hilbert[i]);
	}

	//////////////////////////Calculation of MotionStruct2//////////////////////////

	for (i = 0; i < SysParams->DFTLengthForPSD; i++) {
		Mscan_abs_FFT[i] = (float *) malloc(SysParams->Nbins * sizeof(float));
	}
	for (i = 0; i < SysParams->SpectrogramFreqBins; i++) {
		Pxx2_dB[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
		Pxx2[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
	}
	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert; i++) {

		Pxx2_Hilbert[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
	}

	//PreProcessing
	MacthedFilter2(Mscan2, SysParams);

	//	SlowProcessing2(Mscan2, SysParams); //remove DC
	//	NotchFilter2(Mscan2, SysParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp

	//save the last AvgValue-1 bins of FiftyPrecent for the initial state of the AvgFilter

	for (i = SysParams->AvgValue - 1; i > 0; i--) { //it's done only to Plus and Minus
		Edge2_Plus_2.PrevLastFiftyPrecent[(SysParams->AvgValue - 1) - i] =
				Edge2_Plus_1.FiftyPrecent[SysParams->SpectrogramTimeBins - i];
		Edge2_Minus_2.PrevLastFiftyPrecent[(SysParams->AvgValue - 1) - i] =
				Edge2_Minus_1.FiftyPrecent[SysParams->SpectrogramTimeBins - i];
	}
	//	}
	//Extraction of the curves
	MotionCurveExtraction2(&MotionStruct2, Mscan2, Mscan_abs_FFT, Pxx2_Hilbert,
			Pxx2, Pxx2_dB, &Edge2_2, &Edge2_Plus_2, &Edge2_Minus_2, SysParams);

	///////Create UnitedMotionStruct for the gap interpolation/////////
	UnitedMotionStruct.Edge2 = &Edge2_United;
	UnitedMotionStruct.Edge2_Plus = &Edge2_Plus_United;
	UnitedMotionStruct.Edge2_Minus = &Edge2_Minus_United;
	UnitedMotionStruct.EventStruct = &EventStruct_United;
	UnitedMotionStruct.EventPlusPassEvent = 0;
	UnitedMotionStruct.EventMinusPassEvent = 0;

	//////////////////////Gap Interpolation/////////////////
	GapInterpolation2(MotionStruct0, &MotionStruct1, &MotionStruct2,
			&UnitedMotionStruct, SysParams);

	//////////////////////Feature Extraction/////////////////
	if (SysParams->FirstTimeMotion == 1) {
		SysParams->SpectrogramTimeBins =
				SysParams->SpectrogramTimeBinsTwoMotions;
		FeatureExtractionBasedCurves2(&UnitedMotionStruct, &FeatureSet,
				SysParams);
	} else {
		SysParams->SpectrogramTimeBins =
				SysParams->SpectrogramTimeBinsThreeMotions;
		FeatureExtractionBasedCurves2(&UnitedMotionStruct, &FeatureSet,
				SysParams);

	}

	/////////Classification///////
	memset(MotionDistribution, 0, SysParams->NumOfClasses * sizeof(float));	//initialize with zeros
	if (UnitedMotionStruct.EventPlusPassEvent
			|| UnitedMotionStruct.EventMinusPassEvent) {//there was some event
		y_hat_M = RandomForrestClassifier(&FeatureSet, All_Trees,
				MotionDistribution, RF_Model);
		printf("y_hat %d\n", y_hat_M);
		ClassifierCorrection(&FeatureSet, MotionDistribution, &y_hat_M,
				SVM_Model);

		if ((MotionDistribution[2] < 0.625) && (y_hat_M == 3)) {
			y_hat_M = 1;
			printf("Changed to normal motion from fall, not fully sure\n");
			MotionDistribution[0] = 1;
			MotionDistribution[1] = 0;
			MotionDistribution[2] = 0;
			MotionDistribution[3] = 0;
		}

	} else {	//there is no event, set quasistatic automatically
		y_hat_M = 4;
		MotionDistribution[3] = 1;	//MotionDistribution=[0,0,0,1]
	}

	/////////Save MotionStruct2 into MotionStruc0 for the next loop
	memcpy(MotionStruct0->Edge2->Fmax, MotionStruct2.Edge2->Fmax,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2->FiftyPrecent_Filtered,
			MotionStruct2.Edge2->FiftyPrecent_Filtered,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2->Peak_Filtered,
			MotionStruct2.Edge2->Peak_Filtered,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));

	memcpy(MotionStruct0->Edge2_Plus->Fmax, MotionStruct2.Edge2_Plus->Fmax,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->FiftyPrecent_Filtered,
			MotionStruct2.Edge2_Plus->FiftyPrecent_Filtered,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->Peak_Filtered,
			MotionStruct2.Edge2_Plus->Peak_Filtered,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->SumEnergy_Post,
			MotionStruct2.Edge2_Plus->SumEnergy_Post,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->T1_t, MotionStruct2.Edge2_Plus->T1_t,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));

	memcpy(MotionStruct0->Edge2_Minus->Fmax, MotionStruct2.Edge2_Minus->Fmax,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->FiftyPrecent_Filtered,
			MotionStruct2.Edge2_Minus->FiftyPrecent_Filtered,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->Peak_Filtered,
			MotionStruct2.Edge2_Minus->Peak_Filtered,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->SumEnergy_Post,
			MotionStruct2.Edge2_Minus->SumEnergy_Post,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->T1_t, MotionStruct2.Edge2_Minus->T1_t,
			SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));

	//save the last AvgValue-1 bins of FiftyPrecent for the initial state of the AvgFilter
	for (i = SysParams->AvgValue - 1; i > 0; i--) {	//it's done only to Plus and Minus
		MotionStruct0->Edge2_Plus->PrevLastFiftyPrecent[(SysParams->AvgValue - 1)
														- i] =
																MotionStruct2.Edge2_Plus->FiftyPrecent[SysParams->SpectrogramTimeBinsSingleMotion
																									   - i];
		MotionStruct0->Edge2_Minus->PrevLastFiftyPrecent[(SysParams->AvgValue
				- 1) - i] =
						MotionStruct2.Edge2_Minus->FiftyPrecent[SysParams->SpectrogramTimeBinsSingleMotion
																- i];
	}

	//FREE MEMORY

	for (i = 0; i < SysParams->DFTLengthForPSD; i++) {
		free(Mscan_abs_FFT[i]);
	}
	for (i = 0; i < SysParams->SpectrogramFreqBins; i++) {
		free(Pxx2_dB[i]);
		free(Pxx2[i]);
	}
	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert; i++) {

		free(Pxx2_Hilbert[i]);
	}

	FreeMemory(&Edge2_1, &Edge2_2, &Edge2_Plus_1, &Edge2_Minus_1, &Edge2_Plus_2,
			&Edge2_Minus_2, &UnitedMotionStruct);

	return 0;
}

int GapInterpolation2(Motion_Struct2 *MotionStruct0,
		Motion_Struct2 *MotionStruct1, Motion_Struct2 *MotionStruct2,
		Motion_Struct2 *UnitedMotionStruct, SysParams_Struct* SysParams) {
	int i, FirstGap;
	if (SysParams->FirstTimeMotion == 1) {	//there are only 2 motion structs
		FirstGap = 1;
		//Memory allocation for the new curves with length SpectrogramTimeBinsTwoMotions=2*SpectrogramTimeBins + GapLength
		UnitedMotionStruct->Edge2->Fmax = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Fmax = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Fmax = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));

		UnitedMotionStruct->Edge2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));

		UnitedMotionStruct->Edge2->FiftyPrecent_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered =
				(float *) malloc(
						(SysParams->SpectrogramTimeBinsTwoMotions)
						* sizeof(float));
		UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered =
				(float *) malloc(
						(SysParams->SpectrogramTimeBinsTwoMotions)
						* sizeof(float));

		UnitedMotionStruct->Edge2_Plus->SumEnergy_Post = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->SumEnergy_Post = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));

		UnitedMotionStruct->Edge2_Plus->T1_t = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->T1_t = (float *) malloc(
				(SysParams->SpectrogramTimeBinsTwoMotions) * sizeof(float));

		GapInterpolation_2Curves(MotionStruct1->Edge2->Fmax,
				MotionStruct2->Edge2->Fmax, UnitedMotionStruct->Edge2->Fmax,
				SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->Fmax,
				MotionStruct2->Edge2_Plus->Fmax,
				UnitedMotionStruct->Edge2_Plus->Fmax, SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->Fmax,
				MotionStruct2->Edge2_Minus->Fmax,
				UnitedMotionStruct->Edge2_Minus->Fmax, SysParams, FirstGap);

		GapInterpolation_2Curves(MotionStruct1->Edge2->Peak_Filtered,
				MotionStruct2->Edge2->Peak_Filtered,
				UnitedMotionStruct->Edge2->Peak_Filtered, SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->Peak_Filtered,
				MotionStruct2->Edge2_Plus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Plus->Peak_Filtered, SysParams,
				FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->Peak_Filtered,
				MotionStruct2->Edge2_Minus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Minus->Peak_Filtered, SysParams,
				FirstGap);

		GapInterpolation_2Curves(MotionStruct1->Edge2->FiftyPrecent_Filtered,
				MotionStruct2->Edge2->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered, SysParams,
				FirstGap);
		GapInterpolation_2Curves(
				MotionStruct1->Edge2_Plus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Plus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				SysParams, FirstGap);
		GapInterpolation_2Curves(
				MotionStruct1->Edge2_Minus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Minus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				SysParams, FirstGap);

		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->SumEnergy_Post,
				MotionStruct2->Edge2_Plus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Plus->SumEnergy_Post, SysParams,
				FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->SumEnergy_Post,
				MotionStruct2->Edge2_Minus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post, SysParams,
				FirstGap);

		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->T1_t,
				MotionStruct2->Edge2_Plus->T1_t,
				UnitedMotionStruct->Edge2_Plus->T1_t, SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->T1_t,
				MotionStruct2->Edge2_Minus->T1_t,
				UnitedMotionStruct->Edge2_Minus->T1_t, SysParams, FirstGap);

		SysParams->SpectrogramTimeBins =
				SysParams->SpectrogramTimeBinsTwoMotions;//now the SpectrogramTimebins is longer

		for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {//FiftyPrecent_Filtered corrections, 50precent curve cannot be bigger than the Fmax
			if (UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]
																 > UnitedMotionStruct->Edge2->Fmax[i])
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2->Fmax[i];

			if (UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																	  > UnitedMotionStruct->Edge2_Plus->Fmax[i])
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Plus->Fmax[i];

			if (UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																	   > UnitedMotionStruct->Edge2_Minus->Fmax[i])
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Minus->Fmax[i];
			//			printf("50 two %d %lf\n ", i,
			//					UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]);

		}

	} else {	//there are 3 motion structs
		FirstGap = 1;
		//Memory allocation for the new curves with length SpectrogramTimeBinsThreeMotions=3*SpectrogramTimeBins + 2*GapLength
		UnitedMotionStruct->Edge2->Fmax = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Fmax = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Fmax = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));

		UnitedMotionStruct->Edge2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));

		UnitedMotionStruct->Edge2->FiftyPrecent_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered =
				(float *) malloc(
						(SysParams->SpectrogramTimeBinsThreeMotions)
						* sizeof(float));
		UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered =
				(float *) malloc(
						(SysParams->SpectrogramTimeBinsThreeMotions)
						* sizeof(float));

		UnitedMotionStruct->Edge2_Plus->SumEnergy_Post = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->SumEnergy_Post = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));

		UnitedMotionStruct->Edge2_Plus->T1_t = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->T1_t = (float *) malloc(
				(SysParams->SpectrogramTimeBinsThreeMotions) * sizeof(float));

		//implement gap interpolation between the first two motion structs (MostionStruct0&MostionStruct1)
		GapInterpolation_2Curves(MotionStruct0->Edge2->Fmax,
				MotionStruct1->Edge2->Fmax, UnitedMotionStruct->Edge2->Fmax,
				SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->Fmax,
				MotionStruct1->Edge2_Plus->Fmax,
				UnitedMotionStruct->Edge2_Plus->Fmax, SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->Fmax,
				MotionStruct1->Edge2_Minus->Fmax,
				UnitedMotionStruct->Edge2_Minus->Fmax, SysParams, FirstGap);

		GapInterpolation_2Curves(MotionStruct0->Edge2->Peak_Filtered,
				MotionStruct1->Edge2->Peak_Filtered,
				UnitedMotionStruct->Edge2->Peak_Filtered, SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->Peak_Filtered,
				MotionStruct1->Edge2_Plus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Plus->Peak_Filtered, SysParams,
				FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->Peak_Filtered,
				MotionStruct1->Edge2_Minus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Minus->Peak_Filtered, SysParams,
				FirstGap);

		GapInterpolation_2Curves(MotionStruct0->Edge2->FiftyPrecent_Filtered,
				MotionStruct1->Edge2->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered, SysParams,
				FirstGap);
		GapInterpolation_2Curves(
				MotionStruct0->Edge2_Plus->FiftyPrecent_Filtered,
				MotionStruct1->Edge2_Plus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				SysParams, FirstGap);
		GapInterpolation_2Curves(
				MotionStruct0->Edge2_Minus->FiftyPrecent_Filtered,
				MotionStruct1->Edge2_Minus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				SysParams, FirstGap);

		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->SumEnergy_Post,
				MotionStruct1->Edge2_Plus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Plus->SumEnergy_Post, SysParams,
				FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->SumEnergy_Post,
				MotionStruct1->Edge2_Minus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post, SysParams,
				FirstGap);

		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->T1_t,
				MotionStruct1->Edge2_Plus->T1_t,
				UnitedMotionStruct->Edge2_Plus->T1_t, SysParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->T1_t,
				MotionStruct1->Edge2_Minus->T1_t,
				UnitedMotionStruct->Edge2_Minus->T1_t, SysParams, FirstGap);

		for (i = 0; i < SysParams->SpectrogramTimeBinsTwoMotions; i++) {//FiftyPrecent_Filtered corrections, 50precent curve cannot be bigger than the Fmax
			if (UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]
																 > UnitedMotionStruct->Edge2->Fmax[i])
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2->Fmax[i];

			if (UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																	  > UnitedMotionStruct->Edge2_Plus->Fmax[i])
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Plus->Fmax[i];

			if (UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																	   > UnitedMotionStruct->Edge2_Minus->Fmax[i])
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Minus->Fmax[i];
			//			printf("%d %lf\n",i,UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]);
			//			printf("fmax %d %lf\n",i,UnitedMotionStruct->Edge2_Minus->Fmax[i]);
			//			printf("two UNITED %d %lf\n",i,UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]);
			//			printf("two 50 %d %lf\n", i,
			//					UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]);//HERE CONT

		}

		//now take the result and implement interpolation with the last MostionStruct2
		FirstGap = 0;
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2->Fmax,
				MotionStruct2->Edge2->Fmax, UnitedMotionStruct->Edge2->Fmax,
				SysParams, FirstGap);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->Fmax,
				MotionStruct2->Edge2_Plus->Fmax,
				UnitedMotionStruct->Edge2_Plus->Fmax, SysParams, FirstGap);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Minus->Fmax,
				MotionStruct2->Edge2_Minus->Fmax,
				UnitedMotionStruct->Edge2_Minus->Fmax, SysParams, FirstGap);

		GapInterpolation_2Curves(UnitedMotionStruct->Edge2->Peak_Filtered,
				MotionStruct2->Edge2->Peak_Filtered,
				UnitedMotionStruct->Edge2->Peak_Filtered, SysParams, FirstGap);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->Peak_Filtered,
				MotionStruct2->Edge2_Plus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Plus->Peak_Filtered, SysParams,
				FirstGap);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Minus->Peak_Filtered,
				MotionStruct2->Edge2_Minus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Minus->Peak_Filtered, SysParams,
				FirstGap);

		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered,
				MotionStruct2->Edge2->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered, SysParams,
				FirstGap);
		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Plus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				SysParams, FirstGap);
		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Minus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				SysParams, FirstGap);

		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->SumEnergy_Post,
				MotionStruct2->Edge2_Plus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Plus->SumEnergy_Post, SysParams,
				FirstGap);
		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post,
				MotionStruct2->Edge2_Minus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post, SysParams,
				FirstGap);

		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->T1_t,
				MotionStruct2->Edge2_Plus->T1_t,
				UnitedMotionStruct->Edge2_Plus->T1_t, SysParams, FirstGap);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Minus->T1_t,
				MotionStruct2->Edge2_Minus->T1_t,
				UnitedMotionStruct->Edge2_Minus->T1_t, SysParams, FirstGap);

		for (i = 0; i < SysParams->SpectrogramTimeBinsThreeMotions; i++) {//FiftyPrecent_Filtered corrections, 50precent curve cannot be bigger than the Fmax
			if (UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]
																 > UnitedMotionStruct->Edge2->Fmax[i])
				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2->Fmax[i];

			if (UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																	  > UnitedMotionStruct->Edge2_Plus->Fmax[i])
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Plus->Fmax[i];

			if (UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																	   > UnitedMotionStruct->Edge2_Minus->Fmax[i])
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Minus->Fmax[i];

			//			printf("three peak filtered %d %lf\n", i,
			//					UnitedMotionStruct->Edge2_Plus->Peak_Filtered[i]);//HERE CONT

		}

	}

	return 0;
}

int GapInterpolation_2Curves(float* LeftCurve, float* RightCurve,
		float *IntrpolatedCurve, SysParams_Struct* SysParams, int FirstGap) {
	// we want to find a polynom of 3'rd order y = a3*x^3 + a2*x^2 + a1*x + a0 that best interpolate 2 curves,
	// There are 4 constrains -> and therefore 4 coeff are estimated:
	// 1. p1: continuity in left  curve: the polynom should get the same value in x1 as in curve1(x1):
	// 2. p2: continuity in right curve: the polynom should get the same value in x2 as in curve2(x2):
	//3. p3: same slope/derivative as in curve1 from left:
	// 4. p4: same slope/derivative as in curve2 from right:
	//p1 = a3*x^3 + a2*x^2 + a1*x + a0 @ x1
	// p2 = a3*x^3 + a2*x^2 + a1*x + a0 @ x2
	// p3 = 3a3*x^2 + 2a2*x + a1   + 0  @ x1
	// p4 = 3a3*x^2 + 2a2*x + a1   + 0  @ x2
	int i, m;
	float SumDiffRight = 0, SumDiffLeft = 0;
	float VectorOfValues[4];//0 - LeftVal, 1 - RightVal, 2 - LeftSlope, 3 -RightSlope
	float PolynomCoeffs[4] = { 0 };
	float GapInterpolation[SysParams->GapLength];
	if (FirstGap == 1) {//this is the first interpolation between MotionStruct0 and  MotionStruct1
		VectorOfValues[0] = LeftCurve[SysParams->SpectrogramTimeBinsSingleMotion
									  - 1];
		VectorOfValues[1] = RightCurve[0];

		for (i = 0; i < SysParams->NumSamplesForDerivativeEstimation; i++) {//Calculate mean of the diff
			if (i < SysParams->NumSamplesForDerivativeEstimation - 1)//take 4 elements like in matlab
				SumDiffRight += (RightCurve[i + 1] - RightCurve[i]);
			SumDiffLeft += (LeftCurve[SysParams->SpectrogramTimeBinsSingleMotion
									  - SysParams->NumSamplesForDerivativeEstimation + i]
									  - LeftCurve[SysParams->SpectrogramTimeBinsSingleMotion
												  - SysParams->NumSamplesForDerivativeEstimation + i
												  - 1]);
			//		printf("%d %lf\n",i, LeftCurve[SysParams->SpectrogramTimeBins-SysParams->NumSamplesForDerivativeEstimation+i]);
		}
		VectorOfValues[2] = SumDiffLeft
				/ (SysParams->NumSamplesForDerivativeEstimation);//the diffrenc ein purpose because it's like that in matlab
		VectorOfValues[3] = SumDiffRight
				/ (SysParams->NumSamplesForDerivativeEstimation - 1);

		for (i = 0; i < 4; i++) {	//P=inv(A)*V;
			PolynomCoeffs[i] = SysParams->A12_Inverse[i][0] * VectorOfValues[0]
																			 + SysParams->A12_Inverse[i][1] * VectorOfValues[1]
																															 + SysParams->A12_Inverse[i][2] * VectorOfValues[2]
																																											 + SysParams->A12_Inverse[i][3] * VectorOfValues[3];
		}
		for (m = 0; m < SysParams->GapLength; m++) {//implement polyval  i.e insert the missed timebins in the polynom.
			GapInterpolation[m] = PolynomCoeffs[0]
												* powf((SysParams->SpectrogramTimeBinsSingleMotion + m + 1),
														3)
														+ PolynomCoeffs[1]
																		* powf(
																				(SysParams->SpectrogramTimeBinsSingleMotion
																						+ m + 1), 2)
																						+ PolynomCoeffs[2]
																										* (SysParams->SpectrogramTimeBinsSingleMotion + m
																												+ 1) + PolynomCoeffs[3];

			if (GapInterpolation[m] < SysParams->MinFreq)	//max(gap, minFreq);
				GapInterpolation[m] = SysParams->MinFreq;
		}

		//Insert LeftCurve,GapInterpolation,RightCurve into IntrpolatedCurve
		memcpy(IntrpolatedCurve, LeftCurve,
				SysParams->SpectrogramTimeBinsSingleMotion * sizeof(float));
		memcpy(&IntrpolatedCurve[SysParams->SpectrogramTimeBinsSingleMotion],
				GapInterpolation, SysParams->GapLength * sizeof(float));
		memcpy(
				&IntrpolatedCurve[(SysParams->SpectrogramTimeBinsSingleMotion
						+ SysParams->GapLength)], RightCurve,
						SysParams->SpectrogramTimeBins * sizeof(float));
	} else {//FirstGap=0, this is the second interpolation between MotionStruct1 and  MotionStruct2
		VectorOfValues[0] = LeftCurve[SysParams->SpectrogramTimeBinsTwoMotions
									  - 1];
		VectorOfValues[1] = RightCurve[0];

		for (i = 0; i < SysParams->NumSamplesForDerivativeEstimation; i++) {//Calculate mean of the diff
			if (i < SysParams->NumSamplesForDerivativeEstimation - 1)//take 4 elements like in matlab
				SumDiffRight += (RightCurve[i + 1] - RightCurve[i]);
			SumDiffLeft += (LeftCurve[SysParams->SpectrogramTimeBinsTwoMotions
									  - SysParams->NumSamplesForDerivativeEstimation + i]
									  - LeftCurve[SysParams->SpectrogramTimeBinsTwoMotions
												  - SysParams->NumSamplesForDerivativeEstimation + i
												  - 1]);
		}
		VectorOfValues[2] = SumDiffLeft
				/ (SysParams->NumSamplesForDerivativeEstimation);//the diffrenc ein purpose because it's like that in matlab
		VectorOfValues[3] = SumDiffRight
				/ (SysParams->NumSamplesForDerivativeEstimation - 1);

		for (i = 0; i < 4; i++) {	//P=inv(A)*V;
			PolynomCoeffs[i] = SysParams->A12_Inverse[i][0] * VectorOfValues[0]
																			 + SysParams->A12_Inverse[i][1] * VectorOfValues[1]
																															 + SysParams->A12_Inverse[i][2] * VectorOfValues[2]
																																											 + SysParams->A12_Inverse[i][3] * VectorOfValues[3];
		}
		for (m = 0; m < SysParams->GapLength; m++) {//implement polyval  i.e insert the missed timebins in the polynom.
			GapInterpolation[m] = PolynomCoeffs[0]
												* powf((SysParams->SpectrogramTimeBinsSingleMotion + m + 1),
														3)
														+ PolynomCoeffs[1]
																		* powf(
																				(SysParams->SpectrogramTimeBinsSingleMotion
																						+ m + 1), 2)
																						+ PolynomCoeffs[2]
																										* (SysParams->SpectrogramTimeBinsSingleMotion + m
																												+ 1) + PolynomCoeffs[3];

			if (GapInterpolation[m] < SysParams->MinFreq)	//max(gap, minFreq);
				GapInterpolation[m] = SysParams->MinFreq;
		}

		//Insert LeftCurve,GapInterpolation,RightCurve into IntrpolatedCurve
		memcpy(IntrpolatedCurve, LeftCurve,
				SysParams->SpectrogramTimeBinsTwoMotions * sizeof(float));
		memcpy(&IntrpolatedCurve[SysParams->SpectrogramTimeBinsTwoMotions],
				GapInterpolation, SysParams->GapLength * sizeof(float));
		memcpy(
				&IntrpolatedCurve[(SysParams->SpectrogramTimeBinsTwoMotions
						+ SysParams->GapLength)], RightCurve,
						SysParams->SpectrogramTimeBins * sizeof(float));

	}
	return 0;
}

int ClassifierCorrection(AllFeatures_Struct* FeatureSet,
		float * MotionDistribution, int *y_hat_M, SVM_Struct* SVM_Model) {
	//This function look at specific features in order to correct the classification
	int SNR_thresh = 10;
	if (((FeatureSet->Plus->Max_Black_Curve > 10)
			|| (FeatureSet->Minus->Max_Black_Curve > 10)) && (*y_hat_M == 4)) {
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf(
				"Changed to Normal motion from quasi static, there is a movement\n");
	}
	if ((FeatureSet->Both->SNR_Both < SNR_thresh) && (*y_hat_M != 4)) {
		*y_hat_M = 4;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 1;
		printf("Decision changed to Quasi Static because of low SNR\n");
	}
	if ((FeatureSet->Both->p1 < 0.1) && (*y_hat_M == 3)) {
		*y_hat_M = 1;
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal Motion - probably Standing (P1 < 0.2)\n");
	}
	if (((FeatureSet->Plus->FmaxFpeakMultiplication > 16.5) && (*y_hat_M != 3))) {
		*y_hat_M = 3;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 1;
		MotionDistribution[3] = 0;
		printf(
				"Changed to Fall, high FmaxFpeakMultiplication_Plus ( > 16.5)\n");
	}
	if (((FeatureSet->Plus->yDiff_Fall_Max_FiftyPrecent_Filtered > 6)
			&& (*y_hat_M != 3))) {
		*y_hat_M = 3;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 1;
		MotionDistribution[3] = 0;
		printf("Changed to Fall, high polynomialFeatures_Plus\n");
	}
	if ((FeatureSet->Plus->Edge2_50Precent_PeakToAvg > 2.2)
			&& (*y_hat_M == 2)) {
		*y_hat_M = 3;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 1;
		MotionDistribution[3] = 0;
		printf("Changed to Fall from Sitting,from sitting, big peak2avgplus\n");
	}
	if ((FeatureSet->Plus->Edge2_50Precent_Fmax < 30)
			&& (FeatureSet->Minus->Edge2_50Precent_Fmax < 30)
			&& (*y_hat_M == 3)) {
		*y_hat_M = 1;
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal motion, max blue freq smaller than 30\n");
	}
	if (*y_hat_M == 3)//determine between Fall or False alarms by Garbage(=fast hand movements and get ups form the floor)
		SVMClassifier(FeatureSet, MotionDistribution, y_hat_M, SVM_Model);

	if (MotionDistribution[2] < 0.625 && *y_hat_M == 3) {
		*y_hat_M = 1;
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal motion, not fully sure that it is a fall\n");

	}

	return 0;
}
int SVMClassifier(AllFeatures_Struct* FeatureSet, float * MotionDistribution,
		int *y_hat_M, SVM_Struct *SVM_Model) {
	//This function Normalize the relevant features and calculates the SVM score by
	//score=FeatureSet_normalized*Beta+Bias, if bigger than 0 it's garabge, to avoid misdetect the threshold will be 2
	//The features are: [17	25	26	30	36	41	42]
	float Score = 0;
	int i, NumOfFeatures = 7;
	float NormalizedSVMFeatures[NumOfFeatures];
	NormalizedSVMFeatures[0] =
			(FeatureSet->Plus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq
					- SVM_Model->x_mean[0]) / SVM_Model->x_std[0];
	NormalizedSVMFeatures[1] =
			(FeatureSet->Plus->yDiff_Raise_Max_FiftyPrecent_Filtered
					- SVM_Model->x_mean[1]) / SVM_Model->x_std[1];
	NormalizedSVMFeatures[2] =
			(FeatureSet->Minus->yDiff_Raise_Max_FiftyPrecent_Filtered
					- SVM_Model->x_mean[2]) / SVM_Model->x_std[2];
	NormalizedSVMFeatures[3] =
			(FeatureSet->Minus->yDiff_Raise_Mean_FiftyPrecent_Filtered
					- SVM_Model->x_mean[3]) / SVM_Model->x_std[3];
	NormalizedSVMFeatures[4] = (FeatureSet->Both->p1 - SVM_Model->x_mean[4])
																																											/ SVM_Model->x_std[4];
	NormalizedSVMFeatures[5] = (FeatureSet->Both->maxEventLength
			- SVM_Model->x_mean[5]) / SVM_Model->x_std[5];
	NormalizedSVMFeatures[6] = (FeatureSet->Both->SumEnergy42
			- SVM_Model->x_mean[6]) / SVM_Model->x_std[6];

	for (i = 0; i < NumOfFeatures; i++) {
		Score += SVM_Model->Beta[i] * NormalizedSVMFeatures[i];
	}
	Score += SVM_Model->Bias;	//add the bias

	if (Score > 2 || (Score > 0 && FeatureSet->Both->SumEnergy42 > 0.5)
			|| (Score > 0 && FeatureSet->Plus->Edge2_50Precent_PeakToAvg > 2.83)) {
		*y_hat_M = 1;	//change to normal motion if think that it is garbage
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal motion due to SVMClassifier\n");
	}
	return 0;
}

int RandomForrestClassifier(AllFeatures_Struct* FeatureSet,
		Tree_Struct** All_Trees, float * MotionDistribution,
		RF_Struct* RF_Model) {
	float x[8] = { 0 };	//FeatureSet->NumOFRelevantFeatures=8
	int NaNvalue = -1;
	int y_hat_M = 1;
	int TotalTrees = All_Trees[0]->TotalTrees;	//=8
	int treeIdx, currentPredictor, foundFlag = 0;
	int nextNode, currentNode;	//starting point
	int ProbIdx, NumOfClasses = 4;
	float currentCutPoint, MaxProbabilty;
	//The relevant features in this version is according to MotionClassifierRF_NF_Model8 [11,12,13,18,21,39]
	//The features in MotionClassifierRF_NF_Model9: [11,12,18,21,37,40,41,42], not relevant for this version
	//The decisions are: 1- Normal Motion, 2-Sitting, 3 - Fall, 4- Quasistatic
	//Normalize the features: (feature-mean)/std
	x[0] = (FeatureSet->Plus->Edge2_50Precent_Fmax - RF_Model->x_mean[0])
																																											/ RF_Model->x_std[0];	//#11
	x[1] = (FeatureSet->Minus->Edge2_50Precent_Fmax - RF_Model->x_mean[1])
																																											/ RF_Model->x_std[1];	//#12
	x[2] = (FeatureSet->Plus->Edge2_50Precent_AvgTopFive - RF_Model->x_mean[2])
																																											/ RF_Model->x_std[2];	//#13
	x[3] = (FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq
			- RF_Model->x_mean[3]) / RF_Model->x_std[3];	//#18
	x[4] = (FeatureSet->Both->HilbertRatio - RF_Model->x_mean[4])
																																											/ RF_Model->x_std[4];	//#21
	x[5] = (FeatureSet->Plus->FmaxFpeakMultiplication - RF_Model->x_mean[5])
																																											/ RF_Model->x_std[5];	//#39

	for (treeIdx = 0; treeIdx < TotalTrees; treeIdx++) {
		foundFlag = 0;
		currentNode = 0;
		while (foundFlag == 0) {
			currentPredictor = All_Trees[treeIdx]->CutPredictor[currentNode];//what number of feature
			//			printf("currentPredictor %d\n", currentPredictor);

			if (currentPredictor == NaNvalue) {	//Equivalent to if isnan(currentPredictor)
				foundFlag = 1;	//we are at the end and got decision
				for (ProbIdx = 0; ProbIdx < NumOfClasses; ProbIdx++) {
					*(MotionDistribution + ProbIdx) +=
							All_Trees[treeIdx]->ClassProb[4 * currentNode
														  + ProbIdx];	//sum the final decision probabilities for each tree
					//					printf("%d %d %f\n", treeIdx, ProbIdx,
					//							*(MotionDistribution + ProbIdx));
				}
				break;
			} else {
				currentCutPoint = All_Trees[treeIdx]->CutPoint[currentNode];
				//				printf("cutpoint %lf\n", currentCutPoint);

				if (x[currentPredictor] < currentCutPoint)	//go left
					nextNode = All_Trees[treeIdx]->Children[2 * currentNode];//childern is 2 column array
				else
					//go right
					nextNode =
							All_Trees[treeIdx]->Children[2 * currentNode + 1];

				currentNode = nextNode;
			}
		}

	}

	for (ProbIdx = 0; ProbIdx < NumOfClasses; ProbIdx++) {//calculate mean of the probabilities over all the trees
		*(MotionDistribution + ProbIdx) = *(MotionDistribution + ProbIdx)
																																												/ TotalTrees;
	}
	MaxProbabilty = *(MotionDistribution);
	for (ProbIdx = 1; ProbIdx < NumOfClasses; ProbIdx++) {//find maximum probabilty and y_hat_M
		//		printf("%lf\n",*(MotionDistribution+ProbIdx));
		if (*(MotionDistribution + ProbIdx) > MaxProbabilty) {
			MaxProbabilty=*(MotionDistribution + ProbIdx);
			y_hat_M = ProbIdx + 1;		//final decision
		}
	}

	return y_hat_M;
}

int cmpfunc_for_signal(const void * a, const void * b) {//Comparator function for the quick sort(qsort)
	return (*(float*) b < *(float*) a) - (*(float*) b > *(float*) a);
}

int MotionCurveExtraction2(Motion_Struct2* MotionStruct, float* Mscan[],
		float* Mscan_abs_FFT[], float* Pxx2_Hilbert[], float* Pxx2[],
		float* Pxx2_dB[], Edge2_Struct* Edge2, Edge2_Struct* Edge2_Plus,
		Edge2_Struct* Edge2_Minus, SysParams_Struct* SysParams) {
	int p1,i;
	float T2_dB;
	float* Mscan_PostProcess[SysParams->Nscans];
	for (i=0; i<SysParams->Nscans; i++)
		Mscan_PostProcess[i] = (float *)malloc(SysParams->Nbins * sizeof(float));

	SlowProcessing2(Mscan,Mscan_PostProcess, SysParams); //remove DC
	NotchFilter2(Mscan_PostProcess, SysParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp

	AbsOfFFT2(Mscan_PostProcess, Mscan_abs_FFT, SysParams);
	GET_ROI2(Mscan_abs_FFT, SysParams, &p1);//calculate the region of the bins with max energy
	Spectrogram2(Pxx2, Pxx2_dB, &T2_dB, Mscan_PostProcess, p1, SysParams);//calcualte real spectrogram

	CalcCurves2(Pxx2, Pxx2_dB, &T2_dB, Edge2, SysParams);//CalcRedStarsCurve3 in Matlab


	HilbertSpectrogram2(Pxx2_Hilbert, Pxx2_dB, &T2_dB, Mscan, p1, Edge2_Plus,
			Edge2_Minus, SysParams);//calculate spectrogram after hilbert (complex spectrogram)

	MotionStruct->Edge2 = Edge2;
	MotionStruct->Edge2_Plus = Edge2_Plus;
	MotionStruct->Edge2_Minus = Edge2_Minus;


	for (i = 0; i < SysParams->Nscans; i++) {
		free(Mscan_PostProcess[i]);
	}
	return 0;
}

int FeatureExtractionBasedCurves2(Motion_Struct2* MotionStruct,
		AllFeatures_Struct* FeatureSet, SysParams_Struct* SysParams) {
	int Type;		//1 - ALL ,2 - only Plus, 3- only Minus
	float sum_Edge2_Minus_Fmax = 0, sum_Edge2_Plus_Fmax = 0;
	int i, maxEventLength;

	//Set default values that without event, if there will be event they will be changed
	FeatureSet->Plus->Edge2_50Precent_Fmax = SysParams->MinFreq;//=2.5 Hz   #11;
	FeatureSet->Minus->Edge2_50Precent_Fmax = SysParams->MinFreq; // #12;
	FeatureSet->Plus->Edge2_50Precent_AvgTopFive = SysParams->MinFreq; // #13
	FeatureSet->Plus->Edge2_50Precent_PeakToAvg = 1; //#15
	FeatureSet->Plus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
			SysParams->MinFreq; //#17
	FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
			SysParams->MinFreq; //#18
	FeatureSet->Plus->yDiff_Raise_Max_FiftyPrecent_Filtered = 0; //#25
	FeatureSet->Minus->yDiff_Raise_Max_FiftyPrecent_Filtered = 0; //#26
	FeatureSet->Minus->yDiff_Raise_Mean_FiftyPrecent_Filtered = 0; //#30
	FeatureSet->Plus->Max_Black_Curve = SysParams->MinFreq; //#37
	FeatureSet->Minus->Max_Black_Curve = SysParams->MinFreq; //#38
	FeatureSet->Plus->FmaxFpeakMultiplication = 1; //#39 %SHOULD BE 2.5*2.5/100=0.0625 instead 1 in the origin
	//
	CurveLength2(MotionStruct, MotionStruct->Edge2, SysParams);	//Calculate the event that defined by the peak (black) curve: maxlegnth, start index and end index
	if (MotionStruct->EventStruct->LengthOfMaxEvent
			> SysParams->minEventDuration) {//if the maxlength is bigger than minEventduaration=30
		MotionStruct->EventPlusPassEvent = 1;
		MotionStruct->EventMinusPassEvent = 1;
		Type = 1;	//ALL
		ExtractFeatures2(MotionStruct, FeatureSet, Type, SysParams);
		maxEventLength = MotionStruct->EventStruct->LengthOfMaxEvent;

	} else {
		//try calculate the event based only on Plus
		CurveLength2(MotionStruct, MotionStruct->Edge2_Plus, SysParams);
		if (MotionStruct->EventStruct->LengthOfMaxEvent
				> SysParams->minEventDuration) {//if the maxlength is bigger than minEventduaration=30
			MotionStruct->EventPlusPassEvent = 1;
			Type = 2;	//Only Plus
			ExtractFeatures2(MotionStruct, FeatureSet, Type, SysParams);
			maxEventLength = MotionStruct->EventStruct->LengthOfMaxEvent;
		}
		//try calculate the event based only on Minus
		CurveLength2(MotionStruct, MotionStruct->Edge2_Minus, SysParams);
		if (MotionStruct->EventStruct->LengthOfMaxEvent
				> SysParams->minEventDuration) {//if the maxlength is bigger than minEventduaration=30
			MotionStruct->EventMinusPassEvent = 1;
			Type = 3;	//Only Minus
			ExtractFeatures2(MotionStruct, FeatureSet, Type, SysParams);
			if (MotionStruct->EventStruct->LengthOfMaxEvent > maxEventLength)//for taking the maximum event from plus and minus
				maxEventLength = MotionStruct->EventStruct->LengthOfMaxEvent;
		}

	}

	//#21 HilbertRatio
	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		sum_Edge2_Plus_Fmax += MotionStruct->Edge2_Plus->Fmax[i];
		sum_Edge2_Minus_Fmax += MotionStruct->Edge2_Minus->Fmax[i];
		//		printf("50 filtereed PLUS %d, %f\n", i,
		//				MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]);
	}
	FeatureSet->Both->HilbertRatio = log(
			sum_Edge2_Plus_Fmax / sum_Edge2_Minus_Fmax);
	//

	//#30 SNR Both AND #36  Probability of positive movement
	//	FeatureSet->Both->p1=(MotionStruct->Edge2_Plus->SNR_Linear)/(MotionStruct->Edge2_Plus->SNR_Linear+MotionStruct->Edge2_Minus->SNR_Linear);
	SNRFeature(MotionStruct->Edge2_Plus, MotionStruct->Edge2_Minus, FeatureSet,
			SysParams);
	//

	//#40 maxEventLength
	FeatureSet->Both->maxEventLength = maxEventLength;
	//
	//#42, Energy at high frequencies
	Feature42_2(MotionStruct, FeatureSet);
	//
	return 0;
}

int Feature42_2(Motion_Struct2* MotionStruct, AllFeatures_Struct* FeatureSet) {
	//This function is taking the filtered 50precent curve, zero-padds to 256,perform DFT,normalize and caclulates the energy from min_bin
	//perform to Plus and Minus and take the maximum
	//fast movements like hands have energy at high frequencies so the outcome of this function will be big

	int i, k, FFTLength = 256;
	fftwf_complex out[FFTLength];
	fftwf_plan my_plan;
	int flags = 1;
	float fftin_local_ptr_complex[FFTLength];
	memset(fftin_local_ptr_complex,0,sizeof(fftin_local_ptr_complex));
	float preFFTFiftyPrecentFiltered[FFTLength];
	memset(preFFTFiftyPrecentFiltered, 0, FFTLength * sizeof(float));//initialize with zeros

	int min_bin = 6;//this is the minimum bin that the energy summation will be, equal to 2.7 Hz
	int startIdx = 4;
	float sum_preFFT_FiftyPrecent_Filtered = 0,
			mean_preFFT_FiftyPrecent_Filtered;
	float sum_nonzero_preFFT_FiftyPrecent_Filtered = 0;
	//	_Complex float preFFT_FiftyPrecent_Filtered[FFTLength],
	//	_Complex float  FFT_FiftyPrecent_Filtered[FFTLength / 2 + 1];
	//	_Complex float exp_arg = -I * 2 * M_PI / FFTLength;
	//	memset(preFFT_FiftyPrecent_Filtered, 0, FFTLength * sizeof(_Complex float));//initialize with zeros
	//	float Sum_Energy_Minus, Sum_Energy_Plus, sum_FFT_FiftyPrecent_Filtered;
	float SumEnergyPlus = 0,SumEnergyMinus = 0,sumFFTFiftyPrecentFiltered;
	float absFFTFiftyPrecentFiltered[FFTLength / 2 + 1];
	my_plan = fftwf_plan_dft_r2c_2d(FFTLength,1,fftin_local_ptr_complex,out,flags);



	//	float abs_FFT_FiftyPrecent_Filtered[FFTLength / 2 + 1];
	//Plus
	if (MotionStruct->EventPlusPassEvent == 1) {
		//first zero-pad the signal to Length=256
		if (startIdx < MotionStruct->EventStruct->chosenStart)//max(Event_Plus.chosenStart,5)
			startIdx = MotionStruct->EventStruct->chosenStart;
		for (i = 0; i <= MotionStruct->EventStruct->chosenEnd - startIdx; i++) {
			preFFTFiftyPrecentFiltered[i] =
					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i + startIdx];
			//			preFFT_FiftyPrecent_Filtered[i] =
			//					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i + startIdx];
			sum_preFFT_FiftyPrecent_Filtered +=
					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i + startIdx];
			//			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i+startIdx]!=0)
			//				sum_nonzero_preFFT_FiftyPrecent_Filtered+=1;//sum of non zero elements
		}
		mean_preFFT_FiftyPrecent_Filtered = sum_preFFT_FiftyPrecent_Filtered
				/ FFTLength;
		for (i = 0; i < FFTLength; i++) {//remove DC
			//			preFFT_FiftyPrecent_Filtered[i] -= mean_preFFT_FiftyPrecent_Filtered;
			preFFTFiftyPrecentFiltered[i] = preFFTFiftyPrecentFiltered[i]
																	   - mean_preFFT_FiftyPrecent_Filtered;
			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] != 0)
				sum_nonzero_preFFT_FiftyPrecent_Filtered += 1;//sum of non zero elements
		}
		if (sum_nonzero_preFFT_FiftyPrecent_Filtered != 0) {
			//FFT 256 PTS
			sumFFTFiftyPrecentFiltered = 0;

			//memset(out,0,sizeof(out));
			fftwf_execute_dft_r2c(my_plan,preFFTFiftyPrecentFiltered,out);

			for (k = 0; k <= FFTLength / 2; k++) {	//Length/2 because only one side is wanted
				absFFTFiftyPrecentFiltered[k]=powf(cabs(out[k]),2);
				sumFFTFiftyPrecentFiltered +=
						absFFTFiftyPrecentFiltered[k];
			}
			for (k = min_bin; k <= FFTLength / 2; k++) {
				SumEnergyPlus += absFFTFiftyPrecentFiltered[k]
															/ sumFFTFiftyPrecentFiltered;//normalize and sum from min_bin
			}



			//			sum_FFT_FiftyPrecent_Filtered = 0;
			//
			//			for (k = 0; k <= FFTLength / 2; k++) {	//Length/2 because only one side is wanted
			//				FFT_FiftyPrecent_Filtered[k] = 0;
			//				for (i = 0; i < FFTLength; i++) {
			//					FFT_FiftyPrecent_Filtered[k] += cexp(k * i * exp_arg)
			//																																							* preFFT_FiftyPrecent_Filtered[i];
			//				}
			//				abs_FFT_FiftyPrecent_Filtered[k] = powf(
			//						cabs(FFT_FiftyPrecent_Filtered[k]), 2);
			//				sum_FFT_FiftyPrecent_Filtered +=
			//						abs_FFT_FiftyPrecent_Filtered[k];
			//			}
			//			Sum_Energy_Plus = 0;
			//			for (k = min_bin; k <= FFTLength / 2; k++) {
			//				Sum_Energy_Plus += abs_FFT_FiftyPrecent_Filtered[k]
			//																 / sum_FFT_FiftyPrecent_Filtered;//normalize and sum from min_bin
			//			}
		} else
			SumEnergyPlus = 0;
		//		Sum_Energy_Plus = 0;

	} else
		SumEnergyPlus = 0;
	//		Sum_Energy_Plus = 0;
	//minus
	sum_nonzero_preFFT_FiftyPrecent_Filtered = 0;
	sum_preFFT_FiftyPrecent_Filtered = 0;
	//	memset(preFFT_FiftyPrecent_Filtered, 0, FFTLength * sizeof(_Complex float));//initialize with zeros
	memset(preFFTFiftyPrecentFiltered, 0, FFTLength * sizeof(float));//initialize with zeros
	if (MotionStruct->EventMinusPassEvent == 1) {
		//first zero-pad the signal to Length=256
		if (startIdx < MotionStruct->EventStruct->chosenStart)//max(Event_Plus.chosenStart,5)
			startIdx = MotionStruct->EventStruct->chosenStart;
		for (i = 0; i <= MotionStruct->EventStruct->chosenEnd - startIdx; i++) {
			preFFTFiftyPrecentFiltered[i] =
					MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i + startIdx];
			//			preFFT_FiftyPrecent_Filtered[i] =
			//					MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i
			//																	 + startIdx];
			sum_preFFT_FiftyPrecent_Filtered +=
					MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i
																	 + startIdx];
			//			printf("pre %d %lf\n",i,MotionStruct->Edge2_Minus->FiftyPrecent[i]);
		}
		mean_preFFT_FiftyPrecent_Filtered = sum_preFFT_FiftyPrecent_Filtered
				/ FFTLength;
		for (i = 0; i < FFTLength; i++) {//remove DC
			preFFTFiftyPrecentFiltered[i] -=mean_preFFT_FiftyPrecent_Filtered;

			//			preFFT_FiftyPrecent_Filtered[i] = preFFT_FiftyPrecent_Filtered[i]
			//																		   - mean_preFFT_FiftyPrecent_Filtered;
			if (MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] != 0)
				sum_nonzero_preFFT_FiftyPrecent_Filtered += 1;//sum of non zero elements
		}
		if (sum_nonzero_preFFT_FiftyPrecent_Filtered != 0) {
			//FFT 256 PTS
			sumFFTFiftyPrecentFiltered = 0;

			//memset(out,0,sizeof(out));
			//			my_plan = fftwf_plan_dft_r2c_2d(FFTLength,1,fftin_local_ptr_complex,out,flags);
			fftwf_execute_dft_r2c(my_plan,preFFTFiftyPrecentFiltered,out);

			for (k = 0; k <= FFTLength / 2; k++) {	//Length/2 because only one side is wanted
				absFFTFiftyPrecentFiltered[k]=powf(cabs(out[k]),2);
				sumFFTFiftyPrecentFiltered +=
						absFFTFiftyPrecentFiltered[k];
			}
			for (k = min_bin; k <= FFTLength / 2; k++) {
				SumEnergyMinus += absFFTFiftyPrecentFiltered[k]
															 / sumFFTFiftyPrecentFiltered;//normalize and sum from min_bin
			}



			//
			//
			//			sum_FFT_FiftyPrecent_Filtered = 0;
			//
			//			for (k = 0; k <= FFTLength / 2; k++) {	//Length/2 because only one side is wanted
			//				FFT_FiftyPrecent_Filtered[k] = 0;
			//				for (i = 0; i < FFTLength; i++) {
			//					FFT_FiftyPrecent_Filtered[k] += cexp(k * i * exp_arg)
			//																																							* preFFT_FiftyPrecent_Filtered[i];
			//				}
			//				abs_FFT_FiftyPrecent_Filtered[k] = powf(
			//						cabs(FFT_FiftyPrecent_Filtered[k]), 2);
			//				sum_FFT_FiftyPrecent_Filtered +=
			//						abs_FFT_FiftyPrecent_Filtered[k];
			//			}
			//			Sum_Energy_Minus = 0;
			//			for (k = min_bin; k <= FFTLength / 2; k++) {
			//				Sum_Energy_Minus += abs_FFT_FiftyPrecent_Filtered[k]
			//																  / sum_FFT_FiftyPrecent_Filtered;//normalize and sum from min_bin
			//			}
		} else
			SumEnergyMinus = 0;
		//			Sum_Energy_Minus = 0;
	} else
		SumEnergyMinus = 0;
	//		Sum_Energy_Minus = 0;

	if (SumEnergyPlus > SumEnergyMinus)
		FeatureSet->Both->SumEnergy42 = SumEnergyPlus;
	else
		FeatureSet->Both->SumEnergy42 = SumEnergyMinus;


	fftwf_destroy_plan(my_plan);
	return 0;
}

int ExtractFeatures2(Motion_Struct2* MotionStruct,
		AllFeatures_Struct* FeatureSet, int Type, SysParams_Struct *SysParams) {
	int TimeBin_50Precent_Fmax_Plus = 0, TimeBin_50Precent_Fmax_Minus = 0;
	float Arr_for_sort_Plus[MotionStruct->EventStruct->LengthOfChosenEvent + 1];
	float Edge2_AvgTopFive_Plus = 0;
	float Edge2_noDC_Favg_Plus = 0;
	int num_of_Edge2_noDC_Favg_Plus = 0;
	int i;

	for (i = MotionStruct->EventStruct->chosenStart;
			i <= MotionStruct->EventStruct->chosenEnd; i++) {//run on the timebins in the event

		//		printf("peak filtered %d %lf\n",i ,MotionStruct->Edge2_Plus->Peak_Filtered[i]);

		if (Type == 1 || Type == 2) {			//All or only Plus
			Arr_for_sort_Plus[i - MotionStruct->EventStruct->chosenStart] =
					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];	//for Edge2_AvgTopFive

			//#11 Edge2_50Precent_Fmax_Plus
			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																> FeatureSet->Plus->Edge2_50Precent_Fmax) {
				FeatureSet->Plus->Edge2_50Precent_Fmax =
						MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];
				TimeBin_50Precent_Fmax_Plus = i;			//timebin of Fmax
			}

			//#15 - Edge2+ 50% Peak2Avg: Step 1
			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																> SysParams->MinFreq) {	//calculate mean for fifty_precent_filtered without Fmin bins
				Edge2_noDC_Favg_Plus +=
						MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];
				num_of_Edge2_noDC_Favg_Plus += 1;
			}

			//#37 - max black curve plus
			if (MotionStruct->Edge2_Plus->Peak_Filtered[i]
														> FeatureSet->Plus->Max_Black_Curve)
				FeatureSet->Plus->Max_Black_Curve =
						MotionStruct->Edge2_Plus->Peak_Filtered[i];
			//
		}

		if (Type == 1 || Type == 3) {			//All or only Minus

			//#12 Edge2_50Precent_Fmax_Minus
			if (MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																 > FeatureSet->Minus->Edge2_50Precent_Fmax) {
				FeatureSet->Minus->Edge2_50Precent_Fmax =
						MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i];
				TimeBin_50Precent_Fmax_Minus = i;			//timebin of Fmax

			}
			//			printf("line 1147 %d  ,   %lf\n",i,MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]);

			//

			//#38 - max black curve Minus
			if (MotionStruct->Edge2_Minus->Peak_Filtered[i]
														 > FeatureSet->Minus->Max_Black_Curve)
				FeatureSet->Minus->Max_Black_Curve =
						MotionStruct->Edge2_Minus->Peak_Filtered[i];
			//
		}

	}
	if (Type == 1 || Type == 2) {			//All or only Plus
		//#13 Edge2_50Precent_AvgTopFive_Plus
		qsort(Arr_for_sort_Plus,
				MotionStruct->EventStruct->LengthOfChosenEvent + 1,
				sizeof(float), cmpfunc_for_signal);	//quick sort the 20(=MedianValue) values
		for (i = 0; i < SysParams->topMaxFreq; i++)
			Edge2_AvgTopFive_Plus +=
					Arr_for_sort_Plus[MotionStruct->EventStruct->LengthOfChosenEvent
									  - i];
		FeatureSet->Plus->Edge2_50Precent_AvgTopFive = Edge2_AvgTopFive_Plus
				/ (SysParams->topMaxFreq);
		//		//
		//#15 - Edge2+ 50% Peak2Avg: Step 2
		if (num_of_Edge2_noDC_Favg_Plus == 0)//equivalent to if isempty(Edge2_noDC) in matlab
			num_of_Edge2_noDC_Favg_Plus = SysParams->MinFreq;
		Edge2_noDC_Favg_Plus = Edge2_noDC_Favg_Plus
				/ (num_of_Edge2_noDC_Favg_Plus);
		//
		if (Edge2_noDC_Favg_Plus != 0)	//otherwise Edge2_PeakToAvg_Plus stay 0
			FeatureSet->Plus->Edge2_50Precent_PeakToAvg =
					FeatureSet->Plus->Edge2_50Precent_AvgTopFive
					/ Edge2_noDC_Favg_Plus;
		else
			FeatureSet->Plus->Edge2_50Precent_PeakToAvg = 0;//THINK ABOUT IT AGAIN
		//

		//#17 Edge2_maxPeakFreq_fromTimeBinWithMaxFreq_Plus
		FeatureSet->Plus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
				MotionStruct->Edge2_Plus->Peak_Filtered[TimeBin_50Precent_Fmax_Plus];
		//

		//#25 yDiff_Raise_Max__Edge2_Plus.FiftyPrecent_Filtered & #27 - yDiff_Fall__Max__Edge2_Plus.FiftyPrecent_Filtered
		PolynomialFeatures2(MotionStruct->Edge2_Plus, MotionStruct,
				FeatureSet->Plus);
		//
		//#39 FmaxFpeakMultiplication_Plus
		FeatureSet->Plus->FmaxFpeakMultiplication =
				FeatureSet->Plus->Edge2_50Precent_Fmax
				* (MotionStruct->Edge2_Plus->Peak_Filtered[TimeBin_50Precent_Fmax_Plus])
				* 0.01;
		//
	}
	if (Type == 1 || Type == 3) {		//All or only Minus
		//#18 Edge2_maxPeakFreq_fromTimeBinWithMaxFreq_Minus
		FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
				MotionStruct->Edge2_Minus->Peak_Filtered[TimeBin_50Precent_Fmax_Minus];
		//
		//#26 yDiff_Raise_Max__Edge2_Minus.FiftyPrecent_Filtered, #30 yDiff_Raise_Mean_Edge2_Minus.FiftyPrecent_Filtered
		PolynomialFeatures2(MotionStruct->Edge2_Minus, MotionStruct,
				FeatureSet->Minus);
		//
		//#40 FmaxFpeakMultiplication_Minus
		FeatureSet->Minus->FmaxFpeakMultiplication =
				FeatureSet->Minus->Edge2_50Precent_Fmax
				* (MotionStruct->Edge2_Minus->Peak_Filtered[TimeBin_50Precent_Fmax_Minus])
				* 0.01;
		//
	}

	return 0;
}

int PolynomialFeatures2(Edge2_Struct* Edge2, Motion_Struct2 *MotionStruct,
		Features_Struct *Features) {
	int MedianValue = 20, peakIdx, NumBinsForPolynom;
	int howManySamplesEachSide = 20;
	int timeBin_start = 0;
	int timeBin_stop = MotionStruct->EventStruct->LengthOfChosenEvent;
	int startIdx, stopIdx, i;
	int polynomialOrder = 3; // p4*x^3 + p3*x^2 + p2*x + p1
	double p[polynomialOrder + 1]; //coefficients
	float Edge2_50Precent_MedianFiltered[MotionStruct->EventStruct->LengthOfChosenEvent
										 + MedianValue - 1];
	float MaxValue;
	float yDiff_Raise = 0, yDiff_Fall = 0, yDiff_Raise_mean = 0,
			num_yDiff_Raise_mean = 0;
	//float  yDiff_Fall=0, yDiff_minus_mean=0,num_yDiff_minus_mean=0;//not used

	double *TimeBinsForPolynom, *ValuesForPolynom, yDiff;
	if (howManySamplesEachSide
			> floor(MotionStruct->EventStruct->LengthOfChosenEvent / 2) - 1) //% if the Event is less than 41 samples
		howManySamplesEachSide = floor(
				MotionStruct->EventStruct->LengthOfChosenEvent / 2) - 1;
	MedianFilterFor50Precent(Edge2, MotionStruct,
			Edge2_50Precent_MedianFiltered, MedianValue); //Filter the 50precentfiltered curve
	MaxOfArr(Edge2_50Precent_MedianFiltered, &peakIdx, &MaxValue,
			MotionStruct->EventStruct->LengthOfChosenEvent); //find max idx of 50Precent_MedianFiltered
	//GetMargin in Matlab:
	if (timeBin_start < (peakIdx - howManySamplesEachSide))
		timeBin_start = (peakIdx - howManySamplesEachSide);
	if (timeBin_stop > peakIdx + howManySamplesEachSide)
		timeBin_stop = peakIdx + howManySamplesEachSide;

	if (timeBin_start == 0) {
		startIdx = 0;	//was 1
		stopIdx = 2 * howManySamplesEachSide;	//was+1
	} else if (timeBin_stop == MotionStruct->EventStruct->LengthOfChosenEvent) {
		startIdx = MotionStruct->EventStruct->LengthOfChosenEvent
				- (2 * howManySamplesEachSide) - 1;
		stopIdx = MotionStruct->EventStruct->LengthOfChosenEvent - 1;
	} else {
		startIdx = timeBin_start;
		stopIdx = timeBin_stop;
	}
	//

	NumBinsForPolynom = stopIdx - startIdx + 1;
	TimeBinsForPolynom = (double*) malloc((NumBinsForPolynom) * sizeof(double));//xdata
	ValuesForPolynom = (double*) malloc((NumBinsForPolynom) * sizeof(double));//ydata

	for (i = 0; i < NumBinsForPolynom; i++) {
		TimeBinsForPolynom[i] = startIdx + i + 1;//1 FOR BE LIKE MATLAB????THINKABOUT IT....
		ValuesForPolynom[i] = Edge2_50Precent_MedianFiltered[startIdx + i];
		//		printf("vforpo %d %f\n",i,ValuesForPolynom[i]);
	}
	polyfit(TimeBinsForPolynom, ValuesForPolynom, NumBinsForPolynom,
			polynomialOrder, p);	//MUST TO BE float

	for (i = 0; i <= NumBinsForPolynom; i++) {
		//yDiff_plus = max(yDiff);yDiff_minus = max(-yDiff); yDiff_plus_mean = mean(yDiff(yDiff>=0));yDiff_minus_mean = -mean(yDiff(yDiff<=0));
		if (i < NumBinsForPolynom)
			yDiff = 3 * p[3] * powf(TimeBinsForPolynom[i], 2)
			+ 2 * p[2] * TimeBinsForPolynom[i] + p[1];
		else
			yDiff = 0;

		if (yDiff > yDiff_Raise)	//max(yDiff)
			yDiff_Raise = yDiff;
		if (-yDiff > yDiff_Fall)	//max(-yDiff)
			yDiff_Fall = -yDiff;

		if (yDiff >= 0) {
			yDiff_Raise_mean += yDiff;
			num_yDiff_Raise_mean += 1;
		}

		//		if(-yDiff>yDiff_Fall)
		//			yDiff_Fall=-yDiff;
		//		if (yDiff[i]<=0){//not used
		//			yDiff_minus_mean-=yDiff[i];
		//			num_yDiff_minus_mean+=1;
		//		}
	}
	//#25 yDiff_Raise_Max__Edge2_Plus.FiftyPrecent_Filtered #26yDiff_Raise_Max__Edge2_Minus.FiftyPrecent_Filtered:
	Features->yDiff_Raise_Max_FiftyPrecent_Filtered = yDiff_Raise;
	//
	//#27 yDiff_Fall_Max__Edge2_Plus.FiftyPrecent_Filtered #28yDiff_Fall_Max__Edge2_Minus.FiftyPrecent_Filtered:
	Features->yDiff_Fall_Max_FiftyPrecent_Filtered = yDiff_Fall;
	//
	//#30 - yDiff_Raise_Mean_Edge2_Minus.FiftyPrecent_Filtered
	Features->yDiff_Raise_Mean_FiftyPrecent_Filtered = yDiff_Raise_mean
			/ num_yDiff_Raise_mean;

	//	//#27 - yDiff_Fall__Max__Edge2_Plus.FiftyPrecent_Filtered
	//	Features->yDiff_Fall_Max_FiftyPrecent_Filtered=yDiff_Fall;
	//
	//	Features->yDiff_Fall_Mean_FiftyPrecent_Filtered=yDiff_minus_mean/num_yDiff_minus_mean; not used
	free(TimeBinsForPolynom);
	free(ValuesForPolynom);

	return 0;
}

int SNRFeature(Edge2_Struct* Edge2_Plus, Edge2_Struct* Edge2_Minus,
		AllFeatures_Struct *FeatureSet, SysParams_Struct *SysParams) {
	int i;
	float SumOf_SumEnergy_Post_Plus = 0, MeanOf_SumEnergy_Post_Plus;
	float SumOf_SumEnergy_Post_Minus = 0, MeanOf_SumEnergy_Post_Minus;
	float sum_T1_t_Plus = 0, sum_T1_t_Minus = 0;
	float T1_t_mean_Plus, T1_t_mean_Minus;

	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		SumOf_SumEnergy_Post_Plus += Edge2_Plus->SumEnergy_Post[i];	//for mean calculation
		SumOf_SumEnergy_Post_Minus += Edge2_Minus->SumEnergy_Post[i];
		sum_T1_t_Plus += Edge2_Plus->T1_t[i];	//for mean calculation
		sum_T1_t_Minus += Edge2_Minus->T1_t[i];	//for mean calculation

	}
	T1_t_mean_Plus = (sum_T1_t_Plus / SysParams->SpectrogramTimeBins);
	T1_t_mean_Minus = (sum_T1_t_Minus / SysParams->SpectrogramTimeBins);

	MeanOf_SumEnergy_Post_Plus = SumOf_SumEnergy_Post_Plus
			/ SysParams->SpectrogramTimeBins;
	MeanOf_SumEnergy_Post_Minus = SumOf_SumEnergy_Post_Minus
			/ SysParams->SpectrogramTimeBins;

	Edge2_Plus->SNR_Linear = powf(10,
			(10 * log10(MeanOf_SumEnergy_Post_Plus) - T1_t_mean_Plus) * 0.1);//mean in dB =10*log10(SumOf_SumEnergy_Post_Plus/SpectrogramTimeBins)
	Edge2_Minus->SNR_Linear = powf(10,
			(10 * log10(MeanOf_SumEnergy_Post_Minus) - T1_t_mean_Minus) * 0.1);

	FeatureSet->Both->p1 = (Edge2_Plus->SNR_Linear)
																																											/ (Edge2_Plus->SNR_Linear + Edge2_Minus->SNR_Linear);
	FeatureSet->Both->SNR_Both = (10 * log10(MeanOf_SumEnergy_Post_Plus)
	- T1_t_mean_Plus)
																																											+ (10 * log10(MeanOf_SumEnergy_Post_Minus) - T1_t_mean_Minus);//caluclation for feature SNR_Both

	return 0;
}

int CurveLength2(Motion_Struct2 *MotionStruct, Edge2_Struct* Edge2,
		SysParams_Struct *SysParams) {
	//This function searches for events with minimal length of minEventDuration=30
	//	        | x(n-1)=0    | x(n-1)=1    |
	//	 x(n)=0 | DO NOTHING  | STOP EVENT  |
	//	 x(n)=1 | START EVENT |  DO NOTHING |
	float minValue = 2.512562814070352; //Fmin in Hz
	int i, k, x_before, x_current;
	int size_list = floor(SysParams->SpectrogramTimeBins / 2);
	int startList[size_list];
	int endList[size_list];
	int LengthOfMaxEvent = 0, idxChosenEvent;
	int pointer = 0; //number of the event in startlist & endlist
	float FmaxPerEvent, GlobalFmax = 0;
	int oneHotEdge2_Plus_Peak_Filtered[SysParams->SpectrogramTimeBins + 2]; //+2 for inital state and end state
	memset(oneHotEdge2_Plus_Peak_Filtered, 0,
			(SysParams->SpectrogramTimeBins + 2) * sizeof(int));

	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		if (Edge2->Peak_Filtered[i] > minValue) { //put 1's in all the bins that pass the minValue otherwise stay 0
			oneHotEdge2_Plus_Peak_Filtered[i + 1] = 1;
			//			printf("%d,    %lf\n", i, Edge2->Peak_Filtered[i]);
		} //first and last stayed 0
	}

	for (i = 1; i < SysParams->SpectrogramTimeBins + 2; i++) {
		x_before = oneHotEdge2_Plus_Peak_Filtered[i - 1];
		x_current = oneHotEdge2_Plus_Peak_Filtered[i];

		if ((x_before == 0) && (x_current == 1)) { //start of event
			startList[pointer] = i - 1;
		} else if ((x_before == 1) && (x_current == 0)) { //end of event and
			if (((i - 2) - startList[pointer]) > SysParams->minEventDuration) { //the event is bigger than minduration
				endList[pointer] = i - 2;
				if ((i - 2) - startList[pointer] > LengthOfMaxEvent) {
					LengthOfMaxEvent = (i - 2) - startList[pointer]; //was i-1
					//					idx_of_LengthOfMaxEvent=pointer;
				}
				FmaxPerEvent = 0;
				for (k = startList[pointer]; k <= endList[pointer]; k++) {
					if (Edge2->FiftyPrecent_Filtered[k] > FmaxPerEvent)	//find max of FiftyPrecent_Filtered in the event
						FmaxPerEvent = Edge2->FiftyPrecent_Filtered[k];
				}
				if (FmaxPerEvent > GlobalFmax) {
					GlobalFmax = FmaxPerEvent;
					idxChosenEvent = pointer;
				}
				pointer += 1;
			} else {
				pointer -= 1;//don't take this event in account because it's short
			}
		}
	}
	printf("%d\n", LengthOfMaxEvent);

	if (LengthOfMaxEvent > 0)
		MotionStruct->EventStruct->LengthOfMaxEvent = LengthOfMaxEvent;
	else
		//event not found
		MotionStruct->EventStruct->LengthOfMaxEvent = 0;

	MotionStruct->EventStruct->chosenStart = startList[idxChosenEvent];
	MotionStruct->EventStruct->chosenEnd = endList[idxChosenEvent];
	MotionStruct->EventStruct->LengthOfChosenEvent = endList[idxChosenEvent]
															 - startList[idxChosenEvent];					//was +1

	return 0;
}
int CalcCurvesHilbert2(Pxx2_Plus_Struct *Pxx2_Plus,
		Pxx2_Minus_Struct* Pxx2_Minus, Edge2_Struct* Edge2_Plus,
		Edge2_Struct* Edge2_Minus, SysParams_Struct* SysParams) {//CalcRedStarsCurve3_FirstStep2 in Matlab
	int k, i, idx_FreqBins_Spectogram = 0;
	int MedianType;
	int AllFreqAboveThresh_Plus[SysParams->SpectrogramFreqBinsHilbert / 2],
	AllFreqAboveThresh_Minus[SysParams->SpectrogramFreqBinsHilbert / 2];
	int idx_FreqAboveThresh_Plus, idx_FreqAboveThresh_Minus;
	int MaxPxx2_dB_idx_Plus, MaxPxx2_dB_idx_Minus;
	float freq_increment = 1.25;					//=(Fs/LengthOfDFT);
	float MaxPxx2_dB_Plus, MaxPxx2_dB_Minus;
	//	float SumEnergy_Post_Minus,SumEnergy_Post_Plus;
	//	int Freq_50_PrecentIdx_Plus,Freq_50_PrecentIdx_Minus;
	int flag_50_PrecentIdx_Plus, flag_50_PrecentIdx_Minus;
	float Energy_50_Precent_Plus, Energy_50_Precent_Minus;
	float SpectrogramFreqBinsInHz[SysParams->SpectrogramFreqBinsHilbert / 2];//, Energy_50_Precent;
	for (float i = 0; i < SysParams->SpectrogramFreqBinsHilbert / 2; i++) {	//create the frequency vector inHz
		SpectrogramFreqBinsInHz[idx_FreqBins_Spectogram] = i * freq_increment;
		idx_FreqBins_Spectogram += 1;
	}
	//	memset(Edge2_Plus->Peak,0,174*sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	//	memset(Edge2_Minus->Peak,0,174*sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions

	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		idx_FreqAboveThresh_Plus = 0;
		idx_FreqAboveThresh_Minus = 0;
		MaxPxx2_dB_Plus = Pxx2_Plus->dB[SysParams->Fmin_bin][i];
		MaxPxx2_dB_Minus = Pxx2_Minus->dB[SysParams->Fmin_bin][i];
		MaxPxx2_dB_idx_Plus = SysParams->Fmin_bin;
		MaxPxx2_dB_idx_Minus = SysParams->Fmin_bin;

		for (k = SysParams->Fmin_bin;
				k
				< SysParams->SpectrogramFreqBinsHilbert / 2
				- SysParams->noiseFreqBins; k++) {//frequency range [Fmin_bin,SpectrogramFreqBins-noiseFreqBins]
			if (Pxx2_Plus->dB[k][i] > Edge2_Plus->T1_t[i]) {//if the powfer is passing the threshold
				AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus] = k;
				idx_FreqAboveThresh_Plus += 1;
				//				printf("%d %d\n",idx_FreqAboveThresh_Plus,AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus]);
			}

			if (Pxx2_Plus->dB[k][i] > MaxPxx2_dB_Plus) {//find MaxPxx2_dB_Plus
				MaxPxx2_dB_Plus = Pxx2_Plus->dB[k][i];
				MaxPxx2_dB_idx_Plus = k;
			}

			if (Pxx2_Minus->dB[k][i] > Edge2_Minus->T1_t[i]) {//if the powfer is passing the threshold
				AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus] = k;
				idx_FreqAboveThresh_Minus += 1;
			}

			if (Pxx2_Minus->dB[k][i] > MaxPxx2_dB_Minus) {//find MaxPxx2_dB_Plus
				MaxPxx2_dB_Minus = Pxx2_Minus->dB[k][i];
				MaxPxx2_dB_idx_Minus = k;
			}
		}

		if (idx_FreqAboveThresh_Plus > 0) {	//there were frequencies above the threshold
			Edge2_Plus->maxFreqIdxs[i] =
					AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus - 1];
			Edge2_Plus->maxFreqEnergy[i] =
					Pxx2_Plus->dB[AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus
														  - 1]][i];
			Edge2_Plus->PeakEnergy[i] = MaxPxx2_dB_Plus;
			Edge2_Plus->maxPeakIdxs[i] = MaxPxx2_dB_idx_Plus;
			Edge2_Plus->Fmax[i] =
					SpectrogramFreqBinsInHz[AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus
																	- 1]];

			if (SysParams->truncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Plus->Peak[i + SysParams->MedianValue / 2] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Plus];//first MedianValue/2(=10) bins are 0 for the inital sates of the median filter
			else
				//there is no-zero padding
				Edge2_Plus->Peak[i] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Plus];
			//		Edge2_Plus->SumEnergy[i]=SumEnergy_Plus;
		} else {	//take default values
			Edge2_Plus->maxFreqIdxs[i] = SysParams->Fmin_bin;
			Edge2_Plus->maxFreqEnergy[i] = 0;//should be Pxx2_Plus->dB[Fmin_bin][i]! but in matlab 0..
			Edge2_Plus->PeakEnergy[i] = Pxx2_Plus->dB[SysParams->Fmin_bin][i];
			Edge2_Plus->maxPeakIdxs[i] = SysParams->Fmin_bin;
			Edge2_Plus->Fmax[i] = SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			if (SysParams->truncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Plus->Peak[i + SysParams->MedianValue / 2] =
						SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			else
				//there is no-zero padding
				Edge2_Plus->Peak[i] =
						SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			//		Edge2_Plus->SumEnergy[i]=0;
		}

		if (idx_FreqAboveThresh_Minus > 0) {//there were frequencies above the threshold
			Edge2_Minus->maxFreqIdxs[i] =
					AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus - 1];
			Edge2_Minus->maxFreqEnergy[i] =
					Pxx2_Minus->dB[AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus
															- 1]][i];
			Edge2_Minus->PeakEnergy[i] = MaxPxx2_dB_Minus;
			Edge2_Minus->maxPeakIdxs[i] = MaxPxx2_dB_idx_Minus;
			Edge2_Minus->Fmax[i] =
					SpectrogramFreqBinsInHz[AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus
																	 - 1]];
			if (SysParams->truncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Minus->Peak[i + SysParams->MedianValue / 2] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Minus];//first MedianValue/2(=10) bins are 0 for the inital sates of the median filter
			else
				//there is no-zero padding
				Edge2_Minus->Peak[i] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Minus];//		Edge2_Plus->SumEnergy[i]=SumEnergy_Plus;
		} else {
			Edge2_Minus->maxFreqIdxs[i] = SysParams->Fmin_bin;
			Edge2_Minus->maxFreqEnergy[i] = 0;// should be Pxx2_Minus->dB[Fmin_bin][i];
			Edge2_Minus->PeakEnergy[i] = Pxx2_Minus->dB[SysParams->Fmin_bin][i];
			Edge2_Minus->maxPeakIdxs[i] = SysParams->Fmin_bin;
			Edge2_Minus->Fmax[i] = SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			if (SysParams->truncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Minus->Peak[i + SysParams->MedianValue / 2] =
						SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			else
				//there is no-zero padding
				Edge2_Minus->Peak[i] =
						SpectrogramFreqBinsInHz[SysParams->Fmin_bin];//		Edge2_Plus->SumEnergy[i]=0;
		}

		//		MedianFilter(Edge2_Plus,SpectrogramTimeBins,MedianValue);//median filter on the peak(black) curve with MedianValue=20
		//		MedianFilter(Edge2_Minus,SpectrogramTimeBins,MedianValue);//median filter on the peak(black) curve with MedianValue=20

		if (Edge2_Plus->PeakEnergy[i] > Edge2_Minus->PeakEnergy[i]) {//store the max energy for each timebin
			Edge2_Plus->PeakEnergy_PM[i] = Edge2_Plus->PeakEnergy[i];
			Edge2_Minus->PeakEnergy_PM[i] = Edge2_Plus->PeakEnergy[i];
		} else {
			Edge2_Plus->PeakEnergy_PM[i] = Edge2_Minus->PeakEnergy[i];
			Edge2_Minus->PeakEnergy_PM[i] = Edge2_Minus->PeakEnergy[i];

		}
	}
	MedianType=1;//TRUNCATE
	MedianFilter(Edge2_Plus, SysParams, MedianType);//median filter on the peak(black) curve with MedianValue=20 , MedianType=0 ZERPOAD,MedianType=1 TRUNCATE
	MedianFilter(Edge2_Minus,  SysParams, MedianType);//median filter on the peak(black) curve with MedianValue=20,  MedianType=0 ZERPOAD,MedianType=1 TRUNCATE

	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		Energy_50_Precent_Plus = (Edge2_Plus->PeakEnergy_PM[i]
															+ Edge2_Plus->maxFreqEnergy[i]) / 2;
		Energy_50_Precent_Minus = (Edge2_Minus->PeakEnergy_PM[i]
															  + Edge2_Minus->maxFreqEnergy[i]) / 2;

		//		Freq_50_PrecentIdx_Plus=0;
		//		Freq_50_PrecentIdx_Minus=0;
		flag_50_PrecentIdx_Plus = 0;
		flag_50_PrecentIdx_Minus = 0;

		for (k = Edge2_Plus->maxFreqIdxs[i]; k >= Edge2_Plus->maxPeakIdxs[i];
				k--) {//search backwards the first freq index that passes Energy_50_Precent
			if (Pxx2_Plus->dB[k][i] > Energy_50_Precent_Plus) {
				Edge2_Plus->FiftyPrecent[i] = SpectrogramFreqBinsInHz[k];
				Edge2_Plus->FiftyPrecentIdxs[i] = k;
				flag_50_PrecentIdx_Plus = 1;
				break;
			}
		}
		//		if(Freq_50_PrecentIdx_Plus>Edge2_Plus->maxPeakIdxs[i]){//Freq_50_PrecentIdx_Plus was found
		//			Edge2_Plus->FiftyPrecent[i]=SpectrogramFreqBinsInHz[Freq_50_PrecentIdx_Plus];
		//			Edge2_Plus->FiftyPrecentIdxs[i]=Freq_50_PrecentIdx_Plus;
		//		}
		if (flag_50_PrecentIdx_Plus == 0) {
			Edge2_Plus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2_Plus->Peak_Filtered[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2_Plus->FiftyPrecentIdxs[i] = SysParams->Fmin_bin;
			Edge2_Plus->maxPeakIdxs[i] = SysParams->Fmin_bin;

		}

		for (k = Edge2_Minus->maxFreqIdxs[i]; k >= Edge2_Minus->maxPeakIdxs[i];
				k--) {//search backwards the first freq index that passes Energy_50_Precent
			//			printf("%d %d %lf ",k,i,Pxx2_Minus->dB[k][i]);

			if (Pxx2_Minus->dB[k][i] > Energy_50_Precent_Minus) {
				Edge2_Minus->FiftyPrecent[i] = SpectrogramFreqBinsInHz[k];
				Edge2_Minus->FiftyPrecentIdxs[i] = k;
				flag_50_PrecentIdx_Minus = 1;
				break;
			}
		}
		//		if(Freq_50_PrecentIdx_Minus>Edge2_Minus->maxPeakIdxs[i]){//Freq_50_PrecentIdx_Minus was found
		//			Edge2_Minus->FiftyPrecent[i]=SpectrogramFreqBinsInHz[Freq_50_PrecentIdx_Minus];
		//			Edge2_Minus->FiftyPrecentIdxs[i]=Freq_50_PrecentIdx_Minus;
		//		}
		if (flag_50_PrecentIdx_Minus == 0) {
			Edge2_Minus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2_Minus->Peak_Filtered[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2_Minus->FiftyPrecentIdxs[i] = SysParams->Fmin_bin;
			Edge2_Minus->maxPeakIdxs[i] = SysParams->Fmin_bin;
		}
		//	}

		if (fabs(Edge2_Plus->Peak_Filtered[i])
				<= fabs(SpectrogramFreqBinsInHz[SysParams->Fmin_bin])) {// if the peak has no energy, give the 50% also the same energy
			Edge2_Plus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2_Plus->FiftyPrecentIdxs[i] = SysParams->Fmin_bin;
		}
		if (fabs(Edge2_Minus->Peak_Filtered[i])
				<= fabs(SpectrogramFreqBinsInHz[SysParams->Fmin_bin])) {// if the peak has no energy, give the 50% also the same energy
			Edge2_Minus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2_Minus->FiftyPrecentIdxs[i] = SysParams->Fmin_bin;
		}
		//		SumEnergy_Post_Plus=0;
		//		SumEnergy_Post_Minus=0;
		for (k = SysParams->Fmin_bin; k <= Edge2_Plus->maxPeakIdxs[i]; k++) {
			Edge2_Plus->SumEnergy_Post[i] += Pxx2_Plus->Linear[k][i];//sum the energy until black star curve, in linear
		}
		for (k = SysParams->Fmin_bin; k <= Edge2_Minus->maxPeakIdxs[i]; k++) {
			Edge2_Minus->SumEnergy_Post[i] += Pxx2_Minus->Linear[k][i];	//sum the energy until black star curve, in linear
		}
	}

	AvgFilter2(Edge2_Plus, SysParams->SpectrogramTimeBins, SysParams->AvgValue);//Average filter of the FiftyPrecent(blue) curve with AvgValue=5
	AvgFilter2(Edge2_Minus, SysParams->SpectrogramTimeBins,
			SysParams->AvgValue);//Average filter of the FiftyPrecent(blue) curve with AvgValue=5
	///PRINT HERE FIFTY////////////////
	//	for(i=0;i<SysParams->SpectrogramTimeBins;i++){
	//		printf("###PEAK FILTERED %d %lf\n",i,Edge2_Plus->Peak_Filtered[i] );
	//	}
	return 0;
}

int HilbertSpectrogram2(float* Pxx2_Hilbert[], float* Pxx2_dB[], float *T2_dB,
		float* Mscan[], int p1, Edge2_Struct *Edge2_Plus,
		Edge2_Struct * Edge2_Minus, SysParams_Struct* SysParams) {

	int i, j, m;
	int k;
	int p2 = p1 + SysParams->RangeWide - 1;
	//	_Complex float Signal_for_DFT[SysParams->winLength];

	fftwf_complex *SignalForFFT;
	SignalForFFT= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * SysParams->SpectrogramFreqBinsHilbert);			//pay attention
	memset(SignalForFFT,0,sizeof(fftwf_complex) * SysParams->SpectrogramFreqBinsHilbert);//was in original

	fftwf_complex  *out;
	fftwf_plan p;

	out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * SysParams->SpectrogramFreqBinsHilbert);
	p = fftwf_plan_dft_1d(SysParams->SpectrogramFreqBinsHilbert, SignalForFFT, out, FFTW_FORWARD, FFTW_ESTIMATE);

	int noiseThresh = (5 + 6);
	_Complex float *MscanIQ[SysParams->Nscans];
	//	_Complex float ResultDFT[SysParams->DFTLengthForSpectrogramHilbert];
	//	_Complex float exp_arg = -I * 2 * M_PI
	//			/ SysParams->DFTLengthForSpectrogramHilbert;
	float SumForMeanPlus, SumForMeanMinus;
	//	float* Spectrogram_per_RangeBin[SysParams->SpectrogramFreqBinsHilbert];
	float SpectrogramPerRangeBin[SysParams->SpectrogramFreqBinsHilbert];



	Pxx2_Plus_Struct Pxx2_Plus;
	Pxx2_Minus_Struct Pxx2_Minus;

	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
		Pxx2_Plus.Linear[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));
		Pxx2_Plus.dB[i] = (float *) malloc(
				SysParams->SpectrogramTimeBins * sizeof(float));

		if (i < SysParams->SpectrogramFreqBinsHilbert / 2) {
			Pxx2_Minus.Linear[i] = (float *) malloc(
					SysParams->SpectrogramTimeBins * sizeof(float));
			Pxx2_Minus.dB[i] = (float *) malloc(
					SysParams->SpectrogramTimeBins * sizeof(float));
		}
	}

	//	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert; i++) {
	//		Spectrogram_per_RangeBin[i] = (float *) malloc(
	//				SysParams->SpectrogramTimeBins * sizeof(float));
	//		//				SpectrogramPerRangeBin[i] = (float *) malloc(
	//		//								SysParams->SpectrogramTimeBins * sizeof(float));
	//	}

	for (i = 0; i < SysParams->Nscans; i++) {
		MscanIQ[i] = (_Complex float *) malloc(
				SysParams->Nbins * sizeof(_Complex float));
	}

	Hilbert(Mscan, MscanIQ, SysParams);//get the IQ signal with hilbert transform


	//PreProcess
	SlowProcessingHilbert(MscanIQ,MscanIQ, SysParams); //remove DC
	NotchFilterHilbert(MscanIQ, SysParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp




	for (j = p1; j <= p2; j++) {	//relevant bins [p1,p2]
		for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {//take the next 30 points of slow for the FFT

			for (m = 0; m < SysParams->winLength; m++) {//windowing with hamming
				SignalForFFT[m] = SysParams->Hamming[m]
													 * MscanIQ[m + i * (SysParams->TimeShift)][j];
				//				Signal_for_DFT[m] = SysParams->Hamming[m]
				//													   * MscanIQ[m + i * (SysParams->TimeShift)][j];

			}


			fftwf_execute(p);//shor ttime fft on 30 time bins

			for (k = 0; k < SysParams->SpectrogramFreqBinsHilbert; k++) {

				SpectrogramPerRangeBin[k]= powf(cabsf(out[k]),2);

				if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
					Pxx2_Hilbert[k][i] = 0;
				}
				Pxx2_Hilbert[k][i] += (SpectrogramPerRangeBin[k])
																														/ (SysParams->RangeWide);//Calculate the average for the final spectrogram

				if (k < SysParams->SpectrogramFreqBinsHilbert / 2 + 1) {//fill the Plus part of the spectrogram
					if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
						Pxx2_Plus.Linear[k][i] = 0;
					}
					Pxx2_Plus.Linear[k][i] += (SpectrogramPerRangeBin[k])
																															/ (SysParams->RangeWide);

				}


				if ((k >= SysParams->SpectrogramFreqBinsHilbert / 2)) {	//fill the Minus part of the spectrogram
					if (j == p1) {//if it's the first spectrogram per range bin, put 0 first

						Pxx2_Minus.Linear[SysParams->SpectrogramFreqBinsHilbert - k][i] = 0;

					}
					if (k == SysParams->SpectrogramFreqBinsHilbert / 2) {
						Pxx2_Minus.Linear[0][i] = Pxx2_Plus.Linear[0][i];
					} else {
						Pxx2_Minus.Linear[SysParams->SpectrogramFreqBinsHilbert - k][i] += (SpectrogramPerRangeBin[k])
																																/ (SysParams->RangeWide);//the spectrum stored flipped
					}
				}

			}



			//
			//		for (k = 0; k < SysParams->SpectrogramFreqBinsHilbert; k++) {//Perform short-time fourier transform
			//			ResultDFT[k] = 0;
			//			for (n = 0; n < SysParams->winLength; n++) {
			//				ResultDFT[k] += cexp(k * n * exp_arg) * Signal_for_DFT[n];
			//			}
			//
			//			Spectrogram_per_RangeBin[k][i] = powf(cabs(ResultDFT[k]), 2);
			//			//
			//			if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
			//				Pxx2_Hilbert[k][i] = 0;
			//			}
			//			Pxx2_Hilbert[k][i] += (Spectrogram_per_RangeBin[k][i])
			//																																		/ (SysParams->RangeWide);//Calculate the average for the final spectrogram
			//
			//			if (k < SysParams->SpectrogramFreqBinsHilbert / 2 + 1) {//fill the Plus part of the spectrogram
			//				if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
			//					Pxx2_Plus.Linear[k][i] = 0;
			//				}
			//				Pxx2_Plus.Linear[k][i] += (Spectrogram_per_RangeBin[k][i])
			//																																			/ (SysParams->RangeWide);
			//				//					printf(" p %d %d %lf\n",k,i,Pxx2_Plus.Linear[k][i]);
			//
			//			}
			//			if ((k >= SysParams->SpectrogramFreqBinsHilbert / 2)) {	//fill the Minus part of the spectrogram
			//				if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
			//
			//					Pxx2_Minus.Linear[SysParams->SpectrogramFreqBinsHilbert
			//									  - k][i] = 0;
			//
			//				}
			//				if (k == SysParams->SpectrogramFreqBinsHilbert / 2) {
			//					Pxx2_Minus.Linear[0][i] = Pxx2_Plus.Linear[0][i];
			//				} else {
			//					Pxx2_Minus.Linear[SysParams->SpectrogramFreqBinsHilbert
			//									  - k][i] += (Spectrogram_per_RangeBin[k][i])
			//																																				/ (SysParams->RangeWide);//the spectrum stored flipped
			//				}
			//			}
			//
			//		}


		}
	}

	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		SumForMeanPlus = 0;
		SumForMeanMinus = 0;
		for (k = 0; k < SysParams->SpectrogramFreqBinsHilbert / 2 + 1; k++) {
			Pxx2_Plus.dB[k][i] = 10 * log10(Pxx2_Plus.Linear[k][i]);
			//			printf(" p %d %d %lf %lf %lf\n",i,k,Pxx2_Plus.Linear[k][i],10*log10(Pxx2_Plus.Linear[k][i]),Pxx2_Plus.dB[k][i]);
			if ((k < SysParams->SpectrogramFreqBinsHilbert / 2)) {
				Pxx2_Minus.dB[k][i] = 10 * log10(Pxx2_Minus.Linear[k][i]);
				//				printf(" m %d %d %lf %lf %lf\n",i,k,Pxx2_Minus.Linear[k][i],10*log10(Pxx2_Minus.Linear[k][i]),Pxx2_Minus.dB[k][i]);
			}
			//caclucaltuation of the mean of the last 4 bins for the noise floor estimation
			if (k
					>= (SysParams->SpectrogramFreqBinsHilbert / 2 + 1
							- SysParams->noiseFreqBins)) {
				SumForMeanPlus += Pxx2_Plus.Linear[k][i];
			}
			if (k
					>= (SysParams->SpectrogramFreqBinsHilbert / 2
							- SysParams->noiseFreqBins)) {
				SumForMeanMinus += Pxx2_Minus.Linear[k][i];
			}
		}
		Edge2_Plus->T1_t[i] = 10
				* log10(SumForMeanPlus / SysParams->noiseFreqBins)
		+ noiseThresh;						//WAS PXX2
		Edge2_Minus->T1_t[i] = 10
				* log10(SumForMeanMinus / SysParams->noiseFreqBins)
		+ noiseThresh;

	}

	CalcCurvesHilbert2(&Pxx2_Plus, &Pxx2_Minus, Edge2_Plus, Edge2_Minus,
			SysParams);				//CalcRedStarsCurve3_FirstStep2 in Matlab

	//FREE MEMORY
	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
		free(Pxx2_Plus.Linear[i]);
		free(Pxx2_Plus.dB[i]);

		if (i < SysParams->SpectrogramFreqBinsHilbert / 2) {
			free(Pxx2_Minus.Linear[i]);
			free(Pxx2_Minus.dB[i]);
		}
	}

	for (i = 0; i < SysParams->Nscans; i++) {
		free(MscanIQ[i]);
	}
	//	free(MscanIQ); free it

	//	for (i = 0; i < SysParams->SpectrogramFreqBinsHilbert; i++) {
	//		Spectrogram_per_RangeBin[i] = (float *) malloc(
	//				SysParams->SpectrogramTimeBins * sizeof(float));
	//	}
	//	free(Spectrogram_per_RangeBin); free it


	fftwf_destroy_plan(p);
	fftwf_free(SignalForFFT);
	fftwf_free(out);
	return 0;
}

int Hilbert(float* Mscan[], _Complex float* MscanIQ[], SysParams_Struct* SysParams) {
	//This function is the implemantion of hilbert function in Matlab that based on the article
	//Marple, S. L. �Computing the Discrete-Time Analytic Signal via FFT.�
	int k, n, i;

	float SignalForFFT[SysParams->Nbins];
	fftwf_complex MscanFFT[SysParams->Nbins];
	fftwf_complex MscanIFFT[SysParams->Nbins];

	fftwf_plan my_plan_fft;
	fftwf_plan my_plan_ifft;

	int flags = 1;
	float fftin_local_ptr_complex[SysParams->Nbins];
	memset(fftin_local_ptr_complex,0,sizeof(fftin_local_ptr_complex));

	my_plan_fft = fftwf_plan_dft_r2c_2d(SysParams->Nbins,1,fftin_local_ptr_complex,MscanFFT,flags);

	my_plan_ifft = fftwf_plan_dft_1d(SysParams->Nbins,MscanFFT, MscanIFFT, FFTW_BACKWARD, FFTW_ESTIMATE);

	//	_Complex float Mscan_IFFT[SysParams->Nbins];
	//	memset(Mscan_IFFT,0,sizeof(Mscan_IFFT));

	//	_Complex float exp_arg = -I * 2 * M_PI / SysParams->Nbins;

	for (i = 0; i <  SysParams->Nscans; i++) {//implement FFT, erase the negative frequencies and implement IFFT.

		memcpy(SignalForFFT,Mscan[i],SysParams->Nbins * sizeof(Mscan[0][0]));//prepare the signal for FFT

		fftwf_execute_dft_r2c(my_plan_fft,SignalForFFT,MscanFFT);//FFT
		//		for (n = 0; n < SysParams->Nbins; n++){
		//			printf("%d %d %f + %f\n",i,n,crealf(MscanFFT[n]),cimagf(MscanFFT[n]));
		//
		//		}
		for (k = 1; k < SysParams->Nbins / 2; k++) {
			MscanFFT[k] = 2 * MscanFFT[k]; //in range K=[1,Nscans/2] as the algorithm says and normalize by N

		}
		memset(MscanFFT+(SysParams->Nbins /2+1),0,(SysParams->Nbins / 2)*sizeof(fftwf_complex));//erase the negative part

		fftwf_execute(my_plan_ifft);//IFFT

		for (n = 0; n < SysParams->Nbins; n++){
			MscanIQ[i][n]=conj(MscanIFFT[n])/(SysParams->Nbins);//use memcpy instead
			//printf("%d %d %f + %f\n",i,n,crealf(MscanIFFT[n]),cimagf(MscanIFFT[n]));
		}
		//				memcpy(MscanIQ[i],MscanIFFT/(SysParams->Nbins), sizeof(MscanIFFT));

		//
		//		for (k = 0; k < SysParams->Nbins; k++) {	//FFT
		//			Mscan_FFT[k] = 0;
		//			if (k <= SysParams->Nbins / 2) {//calculate only for the positive frequencies and stay remain 0 the negative
		//				for (n = 0; n < SysParams->Nbins; n++) {//implement one sided FFT in range of k=[0,Nbins/2+1]
		//					Mscan_FFT[k] += cexp(k * n * exp_arg) * Mscan[i][n];
		//
		//				}
		//
		//				if (k >= 1 && k < SysParams->Nbins / 2) {
		//					Mscan_FFT[k] = 2 * Mscan_FFT[k]; //in range K=[1,Nscans/2] as the algorithm says
		//				}
		//			}
		//
		//		}
		//
		//
		//		for (n = 0; n < SysParams->Nbins; n++) { //IFFT
		//			Mscan_IFFT[n] = 0;
		//			for (k = 0; k < SysParams->Nbins; k++) {
		//				Mscan_IFFT[n] += cexp(-k * n * exp_arg) * MscanFFT[k];
		//
		//
		//			}
		////			MscanIQ[i][n] = conj(Mscan_IFFT[n]) / SysParams->Nbins;
		//
		//		}
	}
	fftwf_destroy_plan(my_plan_fft);
	fftwf_destroy_plan(my_plan_ifft);

	return 0;
}

int MedianFilter2(Edge2_Struct *Edge2, int SpectrogramTimeBins, int MedianValue,
		int truncate) { //Filter the black curve (peak curve)
	int i, m;
	//	int delte;
	float ArrForSort[MedianValue];
	float myout[SpectrogramTimeBins];

	if (truncate == 1) { //Computes medians of smaller segments as it reaches the signal edges.
		for (i = 0; i < SpectrogramTimeBins; i++) {

			if (i < (MedianValue / 2)) { //i<10 take the first i+(MedianValue/2) bins
				for (m = 0; m < i + (MedianValue / 2); m++) {
					ArrForSort[m] = Edge2->Peak[m]; //prepare i bins for sorting
				}
				qsort(ArrForSort, i + (MedianValue / 2), sizeof(float),
						cmpfunc_for_signal); //quick sort the i values
				if (i + (MedianValue / 2) % 2 == 0) { //check if even
					Edge2->Peak_Filtered[i] = (*(ArrForSort
							+ (i + MedianValue / 2) / 2 - 1)
							+ *(ArrForSort + (i + MedianValue / 2) / 2)) / 2; //the window is even and therefore the median is the avg

				} else { //i is odd
					Edge2->Peak_Filtered[i] = *(ArrForSort
							+ (i + MedianValue / 2) / 2); //the window length is odd so take the center element

				}
			}

			else if (i <= (SpectrogramTimeBins - MedianValue / 2)) { //i <= (75-10)=65, it means we have 20 bins till the end

				for (m = 0; m < MedianValue; m++) {
					ArrForSort[m] = Edge2->Peak[i + m - MedianValue / 2]; //prepare the 20 bins for sorting
				}

				qsort(ArrForSort, MedianValue, sizeof(float),
						cmpfunc_for_signal); //quick sort the 20(=MedianValue) values
				Edge2->Peak_Filtered[i] = (*(ArrForSort + MedianValue / 2 - 1)
						+ *(ArrForSort + MedianValue / 2)) / 2; //the window is even and therefore the median is the avg
			}
			//			if(i==65)
			//				delte=1;//i> (SpectrogramTimeBins-MedianValue/2
			else if (i > (SpectrogramTimeBins - MedianValue / 2)) {	//i > (75-10)=65, it means we have less than 10 bins till the end, and need to take what is last till the end

				for (m = 0; m < (SpectrogramTimeBins - i + MedianValue / 2);
						m++) {			//m<(75-i-1)
					ArrForSort[m] = Edge2->Peak[i + m - MedianValue / 2];//prepare the SpectrogramTimeBins-i bins for sorting
					//					printf("ARR for sort%d %lf\n", m, ArrForSort[m]);
				}

				qsort(ArrForSort, (SpectrogramTimeBins - i + MedianValue / 2),
						sizeof(float), cmpfunc_for_signal);			//quick sort
				if ((SpectrogramTimeBins - i + MedianValue / 2) % 2 == 0) {	//check if even
					Edge2->Peak_Filtered[i] =
							(*(ArrForSort
									+ (SpectrogramTimeBins - i + MedianValue / 2)
									/ 2 - 1)
									+ *(ArrForSort
											+ (SpectrogramTimeBins - i
													+ MedianValue / 2) / 2))
													/ 2;//the window is even and therefore the median is the avg
					//					printf("even %d %lf\n", i, Edge2->Peak_Filtered[i]);

				} else {			// odd
					Edge2->Peak_Filtered[i] = *(ArrForSort
							+ (SpectrogramTimeBins - i + MedianValue / 2) / 2);	//the window is odd
					//					printf("odd %d %lf\n", i, Edge2->Peak_Filtered[i]);
				}
			}
			//			printf("med %d %lf\n", i, Edge2->Peak_Filtered[i]);

		}


		gsl_vector * input = gsl_vector_alloc(SpectrogramTimeBins);
		gsl_vector * output = gsl_vector_alloc(SpectrogramTimeBins);
		gsl_movstat_workspace * w = gsl_movstat_alloc2(10,9);

		for(i = 0 ; i < SpectrogramTimeBins ; i++)
		{
			gsl_vector_set(input, i, Edge2->Peak[i]);
		}

		gsl_movstat_median(GSL_MOVSTAT_END_TRUNCATE, input, output, w);
		//		myout[0] = gsl_vector_get(output, 0) / 2;
		myout[0] = gsl_vector_get(output, 0);

		for(i = 1 ; i < (SpectrogramTimeBins ) ; i++)
		{
			//				if ((i&0x1) == 0)
			myout[i] =  gsl_vector_get(output, i);
			//				else
			//			tst[i]=gsl_vector_get(output, i);
			//					myout[i] = (gsl_vector_get(output, i) + gsl_vector_get(output, i-1)) / 2;
			printf("value %d : %f\n",i,myout[i] );
		}



	}

	else {			// Considers the signal to be zero beyond the endpoints.
		for (i = 0; i < SpectrogramTimeBins; i++) {
			for (m = 0; m < MedianValue; m++) {
				ArrForSort[m] = Edge2->Peak[i + m];	//prepare the 20 bins for sorting
			}
			qsort(ArrForSort, MedianValue, sizeof(float), cmpfunc_for_signal);//quick sort the 20(=MedianValue) values
			Edge2->Peak_Filtered[i] = (*(ArrForSort + MedianValue / 2 - 1)
					+ *(ArrForSort + MedianValue / 2)) / 2;	//the window is even and therefore the median is the avg
		}


		//		float myout[SpectrogramTimeBins];
		gsl_vector * input = gsl_vector_alloc(SpectrogramTimeBins);
		gsl_vector * output = gsl_vector_alloc(SpectrogramTimeBins);
		//		gsl_movstat_workspace * w = gsl_movstat_alloc(MedianValue+1);
		gsl_movstat_workspace * w = gsl_movstat_alloc2(10,9);

		//int arr[20]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
		for(i = 0 ; i < SpectrogramTimeBins ; i++)
		{
			gsl_vector_set(input, i, Edge2->Peak[i+10]);//need +20
			//				printf("%lf\n", gsl_vector_get(input, j));
		}

		gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);
		//		myout[0] = gsl_vector_get(output, 0) / 2;
		Edge2->Peak_Filtered[0] = gsl_vector_get(output, 0);

		for(i = 1 ; i < (SpectrogramTimeBins ) ; i++)
		{
			//				if ((i&0x1) == 0)
			Edge2->Peak_Filtered[i] =  gsl_vector_get(output, i);
			//				else
			//			tst[i]=gsl_vector_get(output, i);
			//					myout[i] = (gsl_vector_get(output, i) + gsl_vector_get(output, i-1)) / 2;
			//				printf("value %d : %f\n",i,myout[i] );
		}


	}

	//	for(i = 1 ; i < (SpectrogramTimeBins ) ; i++)
	//printf("%d %f\n",i,Edge2->Peak_Filtered[i]);
	return 0;
}

int MedianFilter(Edge2_Struct *Edge2, SysParams_Struct* SysParams,
		int MedianType) { //Filter the black curve (peak curve)
	int i;

	gsl_vector * input = gsl_vector_alloc(SysParams->SpectrogramTimeBins);
	gsl_vector * output = gsl_vector_alloc(SysParams->SpectrogramTimeBins);
	gsl_movstat_workspace * w = gsl_movstat_alloc2(SysParams->MedianValue/2,SysParams->MedianValue/2-1);

	if (MedianType == 1) { //TRUNCATE , Computes medians of smaller segments as it reaches the signal edges.

		for(i = 0 ; i < SysParams->SpectrogramTimeBins ; i++)
		{
			gsl_vector_set(input, i, Edge2->Peak[i]);
		}

		gsl_movstat_median(GSL_MOVSTAT_END_TRUNCATE, input, output, w);

		for(i = 0 ; i < (SysParams->SpectrogramTimeBins ) ; i++)
		{
			Edge2->Peak_Filtered[i] =  gsl_vector_get(output, i);

		}
	}

	else {			// MedianType=0,ZEROPAD, Considers the signal to be zero beyond the endpoints.

		for(i = 0 ; i < SysParams->SpectrogramTimeBins ; i++)
		{
			gsl_vector_set(input, i, Edge2->Peak[i+10]);//need +20
		}

		gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);

		for(i = 0 ; i < (SysParams->SpectrogramTimeBins ) ; i++)
		{
			Edge2->Peak_Filtered[i] =  gsl_vector_get(output, i);

		}
	}
	gsl_movstat_free(w);
	return 0;
}


int MedianFilterFor50Precent(Edge2_Struct *Edge2,
		Motion_Struct2 *MotionStruct, float *Edge2_50Precent_MedianFiltered,
		int MedianValue) {			//Filter the 50precentfiltered curve
	int i;
	int SignalLength = MotionStruct->EventStruct->LengthOfChosenEvent+1;
	gsl_vector * input = gsl_vector_alloc(SignalLength);
	gsl_vector * output = gsl_vector_alloc(SignalLength);
	gsl_movstat_workspace * w = gsl_movstat_alloc2(MedianValue/2,MedianValue/2-1);

	//ZEROPAD, Considers the signal to be zero beyond the endpoints.

	for(i = 0 ; i < SignalLength ; i++)
	{
		gsl_vector_set(input, i, Edge2->FiftyPrecent_Filtered[i + MotionStruct->EventStruct->chosenStart]);//need +20
	}

	gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);

	for(i = 0 ; i < (SignalLength ) ; i++)
	{
		Edge2_50Precent_MedianFiltered[i] =  gsl_vector_get(output, i);
//		printf("%d %f\n",i,Edge2_50Precent_MedianFiltered[i]);
	}

	gsl_movstat_free(w);
	return 0;
}

float MedianFilterFor50Precent2(Edge2_Struct *Edge2,
		Motion_Struct2 *MotionStruct, float *Edge2_50Precent_MedianFiltered,
		int MedianValue) {			//Filter the 50precentfiltered curve
	int i, m;
	int LengthOfSignal = MotionStruct->EventStruct->LengthOfChosenEvent
			+ MedianValue - 1;
	float Arr_for_sort[MedianValue];
	float PreProcess_Arr[LengthOfSignal];//prepare array with MedianValue/2 zeros at beginning and end
	memset(PreProcess_Arr, 0, (LengthOfSignal) * sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions

	for (i = MedianValue / 2;
			i
			<= (MotionStruct->EventStruct->LengthOfChosenEvent
					+ MedianValue / 2); i++) {			//was without =
		PreProcess_Arr[i] = Edge2->FiftyPrecent_Filtered[i - MedianValue / 2
														 + MotionStruct->EventStruct->chosenStart];			////wo -1
	}

	for (i = 0; i <= MotionStruct->EventStruct->LengthOfChosenEvent; i++) {	//was without =
		for (m = 0; m < MedianValue; m++) {
			Arr_for_sort[m] = PreProcess_Arr[i + m];//prepare the 20 bins for sorting
			//			printf("arr dor sort  50 %d   %f\n",i,Arr_for_sort[m]);
		}
		qsort(Arr_for_sort, MedianValue, sizeof(float), cmpfunc_for_signal);//quick sort the 20(=MedianValue) values
		Edge2_50Precent_MedianFiltered[i] = (*(Arr_for_sort + MedianValue / 2
				- 1) + *(Arr_for_sort + MedianValue / 2)) / 2;//the window is even and therefore the median is the avg
		//	printf("med 50 %d   %f\n",i,Edge2_50Precent_MedianFiltered[i]);
	}


	int SignalLength = MotionStruct->EventStruct->LengthOfChosenEvent+1;


	gsl_vector * input = gsl_vector_alloc(SignalLength);
	gsl_vector * output = gsl_vector_alloc(SignalLength);
	gsl_movstat_workspace * w = gsl_movstat_alloc2(MedianValue/2,MedianValue/2-1);

	//ZEROPAD, Considers the signal to be zero beyond the endpoints.

	for(i = 0 ; i < SignalLength ; i++)
	{
		gsl_vector_set(input, i, Edge2->FiftyPrecent_Filtered[i + MotionStruct->EventStruct->chosenStart]);//need +20
	}

	gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);

	for(i = 0 ; i < (SignalLength ) ; i++)
	{
		Edge2_50Precent_MedianFiltered[i] =  gsl_vector_get(output, i);
		printf("%d %f\n",i,Edge2_50Precent_MedianFiltered[i]);
	}


	return 0;
}

int SlowProcessing2(float* Mscan[],float* Mscan_PostProcess[], SysParams_Struct* SysParams) {
	int i, j, WindowLength = 10;
	float Alpha = 0.1;
	float AdaptiveDC[SysParams->Nscans + WindowLength];
	float MeanOf10FirstTimeBins;
	for (j = 0; j < SysParams->Nbins; j++) {
		MeanOf10FirstTimeBins = 0;
		for (int k = 0; k < WindowLength; k++) {//Calculate mean of 10 first bins for the initial conditions
			MeanOf10FirstTimeBins += Mscan[k][j];
		}
		MeanOf10FirstTimeBins = MeanOf10FirstTimeBins / WindowLength;
		AdaptiveDC[0] = (1 - Alpha) * MeanOf10FirstTimeBins
				+ Alpha * Mscan[0][j];

		for (i = 1; i < SysParams->Nscans + WindowLength; i++) {
			if (i < SysParams->Nscans - 1) {
				AdaptiveDC[i] = (1 - Alpha) * AdaptiveDC[i - 1]
														 + Alpha * Mscan[i][j];
			} else {		//if in the last 10 time bins don't change the value
				AdaptiveDC[i] = AdaptiveDC[i - 1];
			}
		}
		for (i = 0; i < SysParams->Nscans; i++) {
			//			Mscan[i][j] -= AdaptiveDC[i + WindowLength];
			Mscan_PostProcess[i][j]=Mscan[i][j] - AdaptiveDC[i + WindowLength];
			//						printf("%d %d %lf\n",i, j,Mscan[i][j]);

		}
	}

	return 0;
}


int SlowProcessingHilbert(_Complex float* Mscan[],_Complex float* Mscan_PostProcess[], SysParams_Struct* SysParams) {
	int i, j, WindowLength = 10;
	_Complex float Alpha = 0.1;
	_Complex float AdaptiveDC[SysParams->Nscans + WindowLength];
	_Complex float MeanOf10FirstTimeBins;
	for (j = 0; j < SysParams->Nbins; j++) {
		MeanOf10FirstTimeBins = 0;
		for (int k = 0; k < WindowLength; k++) {//Calculate mean of 10 first bins for the initial conditions
			MeanOf10FirstTimeBins += Mscan[k][j];
		}
		MeanOf10FirstTimeBins = MeanOf10FirstTimeBins / WindowLength;
		AdaptiveDC[0] = (1 - Alpha) * MeanOf10FirstTimeBins
				+ Alpha * Mscan[0][j];

		for (i = 1; i < SysParams->Nscans + WindowLength; i++) {
			if (i < SysParams->Nscans - 1) {
				AdaptiveDC[i] = (1 - Alpha) * AdaptiveDC[i - 1]
														 + Alpha * Mscan[i][j];
			} else {		//if in the last 10 time bins don't change the value
				AdaptiveDC[i] = AdaptiveDC[i - 1];
			}
		}
		for (i = 0; i < SysParams->Nscans; i++) {
			//			Mscan[i][j] -= AdaptiveDC[i + WindowLength];
			Mscan_PostProcess[i][j]=Mscan[i][j] - AdaptiveDC[i + WindowLength];
			//						printf("%d %d %lf\n",i, j,Mscan[i][j]);

		}
	}

	return 0;
}

int NotchFilter2(float* Mscan[], SysParams_Struct* SysParams) {

	//first num_concatenate_bins the first 50 bins to the beginning of the Mscan for the initial conditions
	//then filter the slows with notch filter at 50 and 100 Hz.
	int i, j, num_concatenate_scans = 50;
	float concatenated_slow[SysParams->Nscans + num_concatenate_scans],
	filter_result_50[SysParams->Nscans + num_concatenate_scans],
	filter_result_100[SysParams->Nscans + num_concatenate_scans];
	float a_50[3], a_100[3], b_50[3], b_100[3];
	//filter coefficients:
	a_50[0] = 1;
	a_50[1] = -0.546618432826014;
	a_50[2] = 0.768894406379384;
	b_50[0] = 0.884447203189692;
	b_50[1] = -0.546618432826014;
	b_50[2] = 0.884447203189692;
	a_100[0] = 1;
	a_100[1] = 1.43106563601571;
	a_100[2] = 0.768894406379384;
	b_100[0] = 0.884447203189692;
	b_100[1] = 1.43106563601571;
	b_100[2] = 0.884447203189692;

	for (j = 0; j < SysParams->Nbins; j++) {	//run on the columns
		for (i = 0; i < SysParams->Nscans + num_concatenate_scans; i++) {
			if (i < num_concatenate_scans) {//concatenate the first 50 nscans
				concatenated_slow[i] = Mscan[i][j];
			} else {
				concatenated_slow[i] = Mscan[i - num_concatenate_scans][j];
			}
			//notch filter @50Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_50[0] = b_50[0] * concatenated_slow[0] / a_50[0];
			} else if (i == 1) {
				filter_result_50[1] = (b_50[0] * concatenated_slow[1]
																   + b_50[1] * concatenated_slow[0]
																								 - a_50[1] * filter_result_50[0]) / a_50[0];
			} else {

				filter_result_50[i] = (b_50[0] * concatenated_slow[i]
																   + b_50[1] * concatenated_slow[i - 1]
																								 + b_50[2] * concatenated_slow[i - 2]
																															   - a_50[1] * filter_result_50[i - 1]
																																							- a_50[2] * filter_result_50[i - 2]) / a_50[0];
			}
			//notch filter @100Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_100[0] = b_100[0] * filter_result_50[0]
																   / a_100[0];
			} else if (i == 1) {
				filter_result_100[1] = (b_100[0] * filter_result_50[1]
																	+ b_100[1] * filter_result_50[0]
																								  - a_100[1] * filter_result_100[0]) / a_100[0];
			} else {

				filter_result_100[i] = (b_100[0] * filter_result_50[i]
																	+ b_100[1] * filter_result_50[i - 1]
																								  + b_100[2] * filter_result_50[i - 2]
																																- a_100[1] * filter_result_100[i - 1]
																																							   - a_100[2] * filter_result_100[i - 2]) / a_100[0];
			}
			//			printf("%d %d %lf\n",i, j,filter_result_100[i]);
			if (i >= num_concatenate_scans) {	//crop the the relevant Nscans
				Mscan[i - num_concatenate_scans][j] = filter_result_100[i];
				//								printf("%d %d %lf \n",i-num_concatenate_scans, j,Mscan[i-num_concatenate_scans][j]);
			}
		}
	}

	return 0;
}


int NotchFilterHilbert(_Complex float* Mscan[], SysParams_Struct* SysParams) {

	//first num_concatenate_bins the first 50 bins to the beginning of the Mscan for the initial conditions
	//then filter the slows with notch filter at 50 and 100 Hz.
	int i, j, num_concatenate_scans = 50;
	_Complex float concatenated_slow[SysParams->Nscans + num_concatenate_scans],
	filter_result_50[SysParams->Nscans + num_concatenate_scans],
	filter_result_100[SysParams->Nscans + num_concatenate_scans];
	_Complex float a_50[3], a_100[3], b_50[3], b_100[3];
	//filter coefficients:
	a_50[0] = 1;
	a_50[1] = -0.546618432826014;
	a_50[2] = 0.768894406379384;
	b_50[0] = 0.884447203189692;
	b_50[1] = -0.546618432826014;
	b_50[2] = 0.884447203189692;
	a_100[0] = 1;
	a_100[1] = 1.43106563601571;
	a_100[2] = 0.768894406379384;
	b_100[0] = 0.884447203189692;
	b_100[1] = 1.43106563601571;
	b_100[2] = 0.884447203189692;

	for (j = 0; j < SysParams->Nbins; j++) {	//run on the columns
		for (i = 0; i < SysParams->Nscans + num_concatenate_scans; i++) {
			if (i < num_concatenate_scans) {//concatenate the first 50 nscans
				concatenated_slow[i] = Mscan[i][j];
			} else {
				concatenated_slow[i] = Mscan[i - num_concatenate_scans][j];
			}
			//notch filter @50Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_50[0] = b_50[0] * concatenated_slow[0] / a_50[0];
			} else if (i == 1) {
				filter_result_50[1] = (b_50[0] * concatenated_slow[1]
																   + b_50[1] * concatenated_slow[0]
																								 - a_50[1] * filter_result_50[0]) / a_50[0];
			} else {

				filter_result_50[i] = (b_50[0] * concatenated_slow[i]
																   + b_50[1] * concatenated_slow[i - 1]
																								 + b_50[2] * concatenated_slow[i - 2]
																															   - a_50[1] * filter_result_50[i - 1]
																																							- a_50[2] * filter_result_50[i - 2]) / a_50[0];
			}
			//notch filter @100Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_100[0] = b_100[0] * filter_result_50[0]
																   / a_100[0];
			} else if (i == 1) {
				filter_result_100[1] = (b_100[0] * filter_result_50[1]
																	+ b_100[1] * filter_result_50[0]
																								  - a_100[1] * filter_result_100[0]) / a_100[0];
			} else {

				filter_result_100[i] = (b_100[0] * filter_result_50[i]
																	+ b_100[1] * filter_result_50[i - 1]
																								  + b_100[2] * filter_result_50[i - 2]
																																- a_100[1] * filter_result_100[i - 1]
																																							   - a_100[2] * filter_result_100[i - 2]) / a_100[0];
			}
			//			printf("%d %d %lf\n",i, j,filter_result_100[i]);
			if (i >= num_concatenate_scans) {	//crop the the relevant Nscans
				Mscan[i - num_concatenate_scans][j] = filter_result_100[i];
				//								printf("%d %d %lf \n",i-num_concatenate_scans, j,Mscan[i-num_concatenate_scans][j]);
			}
		}
	}

	return 0;
}

int AbsOfFFT2(float* Mscan[], float* Mscan_abs_FFT[],
		SysParams_Struct* SysParams) {
	int k, n, j;
	_Complex float Mscan_FFT[SysParams->DFTLengthForPSD];

	fftwf_complex MscanFFT[SysParams->DFTLengthForPSD];

	fftwf_plan my_plan;
	int flags = 1;
	float SignalForFFT[SysParams->Nscans];
	float fftin_local_ptr_complex[SysParams->DFTLengthForPSD];
	memset(fftin_local_ptr_complex,0,sizeof(fftin_local_ptr_complex));
	my_plan = fftwf_plan_dft_r2c_2d(SysParams->Nscans,1,fftin_local_ptr_complex,MscanFFT,flags);


	for (j = 0; j < SysParams->Nbins; j++) {
		for (n = 0; n < SysParams->Nscans; n++){
			SignalForFFT[n]=Mscan[n][j];
		}
		fftwf_execute_dft_r2c(my_plan,SignalForFFT,MscanFFT);

		for (k = 0; k <SysParams->DFTLengthForPSD; k++) {	//Length/2 because only one side is wanted
			Mscan_abs_FFT[k][j]=powf(cabs(MscanFFT[k]),2);
		}


		////////////////////////
		for (k = 0; k < SysParams->DFTLengthForPSD; k++) {

			Mscan_FFT[k] = 0;
			for (n = 0; n < SysParams->Nscans; n++) {//implement one sided FFT in range of k=[0,Nscans/2+1]
				Mscan_FFT[k] += cexp(-I * k * n * 2 * M_PI / SysParams->Nscans)
																																														* Mscan[n][j];
			}
			Mscan_abs_FFT[k][j] = powf(cabs(Mscan_FFT[k]), 2);
		}
	}
	fftwf_destroy_plan(my_plan);
	return 0;
}

int MacthedFilter(float* Mscan[], int Nscans, int Nbins) {
	int Filter_Length = 19;
	float filter_coeffs[19];
	float InverseFilterGain = 0.599264554413850;
	float filtered_result[Nbins];
	int i, j, running_idx;

	filter_coeffs[0] = 0.000408;
	filter_coeffs[1] = -0.041049471281062;
	filter_coeffs[2] = 0.008649261539219;
	filter_coeffs[3] = 0.229369;
	filter_coeffs[4] = -0.0563189130046412;
	filter_coeffs[5] = -0.520325728953466;
	filter_coeffs[6] = 0.165011;
	filter_coeffs[7] = 0.892380443834730;
	filter_coeffs[8] = -0.067669198947082;
	filter_coeffs[9] = -0.9717;
	filter_coeffs[10] = 0.156922538937629;
	filter_coeffs[11] = 0.717920067065040;
	filter_coeffs[12] = -0.128898;
	filter_coeffs[13] = -0.349983210110945;
	filter_coeffs[14] = 0.0915642458163692;
	filter_coeffs[15] = 0.090303;
	filter_coeffs[16] = -0.022514346739832;
	filter_coeffs[17] = -0.013823220906457;
	filter_coeffs[18] = 0.003896;

	for (i = 0; i < Nscans; i++) {
		for (j = 0; j < Nbins; j++) {

			filtered_result[j] = 0;
			if (j < Filter_Length) {//the case of filtration of the first 19 bins
				for (running_idx = 0; running_idx <= j; running_idx++)
					filtered_result[j] += Mscan[i][j - running_idx]
												   * filter_coeffs[running_idx];

			} else {
				for (running_idx = 0; running_idx < Filter_Length;
						running_idx++)
					filtered_result[j] += Mscan[i][j - running_idx]
												   * filter_coeffs[running_idx];
			}

		}
		for (j = 0; j < Nbins; j++) {
			Mscan[i][j] = filtered_result[j] * InverseFilterGain;
			;
			//			printf("%d %d %lf\n",i, j,Mscan[i][j]);
		}
	}

	return 0;
}

int MacthedFilter2(float* Mscan[], SysParams_Struct* SysParams) {
	int Filter_Length = 19;
	float filter_coeffs[19];
	float InverseFilterGain = 0.599264554413850;
	float filtered_result[SysParams->Nbins];
	int i, j, running_idx;

	filter_coeffs[0] = 0.000408;
	filter_coeffs[1] = -0.041049471281062;
	filter_coeffs[2] = 0.008649261539219;
	filter_coeffs[3] = 0.229369;
	filter_coeffs[4] = -0.0563189130046412;
	filter_coeffs[5] = -0.520325728953466;
	filter_coeffs[6] = 0.165011;
	filter_coeffs[7] = 0.892380443834730;
	filter_coeffs[8] = -0.067669198947082;
	filter_coeffs[9] = -0.9717;
	filter_coeffs[10] = 0.156922538937629;
	filter_coeffs[11] = 0.717920067065040;
	filter_coeffs[12] = -0.128898;
	filter_coeffs[13] = -0.349983210110945;
	filter_coeffs[14] = 0.0915642458163692;
	filter_coeffs[15] = 0.090303;
	filter_coeffs[16] = -0.022514346739832;
	filter_coeffs[17] = -0.013823220906457;
	filter_coeffs[18] = 0.003896;

	for (i = 0; i < SysParams->Nscans; i++) {
		for (j = 0; j < SysParams->Nbins; j++) {

			filtered_result[j] = 0;
			if (j < Filter_Length) {//the case of filtration of the first 19 bins
				for (running_idx = 0; running_idx <= j; running_idx++)
					filtered_result[j] += Mscan[i][j - running_idx]
												   * filter_coeffs[running_idx];

			} else {
				for (running_idx = 0; running_idx < Filter_Length;
						running_idx++)
					filtered_result[j] += Mscan[i][j - running_idx]
												   * filter_coeffs[running_idx];
			}

		}
		for (j = 0; j < SysParams->Nbins; j++) {
			Mscan[i][j] = filtered_result[j] * InverseFilterGain;
			;
			//			printf("%d %d %lf\n",i, j,Mscan[i][j]);
		}
	}

	return 0;
}

int GET_ROI2(float *Mscan_abs_FFT[], SysParams_Struct* SysParams, int *p1) {
	int FrequncyRange = ceil(SysParams->DFTLengthForPSD / 2);//Rangepowfer is calucalted over this range of frequencies
	float Rangepowfer[SysParams->Nbins];
	float y[SysParams->Nbins - SysParams->RangeWide + 1], max_y;
	int i, j, m;

	for (j = 0; j < SysParams->Nbins; j++) {
		Rangepowfer[j] = 0;
		for (i = 0; i <= FrequncyRange; i++) {
			Rangepowfer[j] += Mscan_abs_FFT[i][j];
		}
	}

	for (j = 0; j < SysParams->Nbins - SysParams->RangeWide + 1; j++) {
		y[j] = 0;
		for (m = j; m <= SysParams->RangeWide - 1 + j; m++) {
			y[j] += Rangepowfer[m];
		}
		if (j == 0) {
			max_y = y[0];
			*p1 = 0;
		}
		if (y[j] > max_y) {
			max_y = y[j];
			*p1 = j;	//p1=idx of max_y
		}
	}
	return 0;
}

int Spectrogram2(float* Pxx2[], float* Pxx2_dB[], float *T2_dB, float* Mscan[],


		int p1, SysParams_Struct* SysParams) {


	fftwf_complex out[SysParams->DFTLengthForSpectrogram];

	float SignalForFFT[SysParams->DFTLengthForSpectrogram];
	memset(SignalForFFT,0,sizeof(SignalForFFT));//was in original

	float SpectrogramPerRangeBin[SysParams->SpectrogramFreqBins][SysParams->SpectrogramTimeBins];


	fftw_plan my_plan;
	int flags = 1;
	float fftin_local_ptr_complex[SysParams->DFTLengthForSpectrogram];
	memset(fftin_local_ptr_complex,0,sizeof(fftin_local_ptr_complex));

	int i, j, m;
	int k, n, p2;
	float Signal_for_DFT[SysParams->winLength];

	int SizeForMean = SysParams->SpectrogramTimeBins
			* (SysParams->noiseFreqBins);
	//	float Spectrogram_per_RangeBin[SysParams->DFTLengthForSpectrogram / 2 + 1][SysParams->SpectrogramTimeBins];
	float T2_Linear = 0;	//average of all the spectrogram in linear
	//	_Complex float ResultDFT[SysParams->DFTLengthForSpectrogram];
	//	_Complex float exp_arg = -I * 2 * M_PI / SysParams->DFTLengthForSpectrogram;

	p2 = p1 + SysParams->RangeWide - 1;
	for (j = p1; j <= p2; j++) {	//relevant bins [p1,p2]
		for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {//take the next 30 points of slow for the FFT
			for (m = 0; m < SysParams->winLength; m++) {//windowing with hamming
				//				Signal_for_DFT[m] = SysParams->Hamming[m]
				//													   * Mscan[m + i * (SysParams->TimeShift)][j];
				SignalForFFT[m] = SysParams->Hamming[m]
													 * Mscan[m + i * (SysParams->TimeShift)][j];
			}

			//FFT
			//	memset(out,0,sizeof(out));//was in original
			my_plan = fftwf_plan_dft_r2c_2d(SysParams->DFTLengthForSpectrogram,1,fftin_local_ptr_complex,out,flags);
			fftwf_execute_dft_r2c(my_plan,SignalForFFT,out);
			for (k = 0; k < SysParams->SpectrogramFreqBins; k++) {//Perform short-time fourier transform take in range [0,pi]

				SpectrogramPerRangeBin[k][i] = powf(cabs(out[k]),2);
				if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
					Pxx2[k][i] = 0;
				}
				Pxx2[k][i] += (SpectrogramPerRangeBin[k][i])
																																																	/ (SysParams->RangeWide);//Calculate the average for the final spectrogram
			}


			//			for (k = 0; k < SysParams->SpectrogramFreqBins; k++) {//Perform short-time fourier transform take in range [0,pi]
			//				ResultDFT[k] = 0;
			//				for (n = 0; n < SysParams->winLength; n++) {
			//					ResultDFT[k] += cexp(n * k * exp_arg) * Signal_for_DFT[n];
			//				}
			//				Spectrogram_per_RangeBin[k][i] = powf(cabs(ResultDFT[k]), 2);
			//				//							printf(" %d %d %lf\n",k,i,Spectrogram_per_RangeBin[k][i]);
			//
			//				if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
			//					Pxx3[k][i] = 0;
			//				}
			//				Pxx3[k][i] += (Spectrogram_per_RangeBin[k][i])
			//															/ (SysParams->RangeWide);//Calculate the average for the final spectrogram
			//				//							printf("%d %d    %f\n",k,i,Pxx2[k][i]);
			//			}
		}
	}
	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {
		//		for(k=SpectrogramFreqBins-noiseFreqBins;k<SpectrogramFreqBins;k++){
		for (k = 0; k < SysParams->SpectrogramFreqBins; k++) {

			if ((k >= SysParams->SpectrogramFreqBins - SysParams->noiseFreqBins)
					&& (k < SysParams->SpectrogramFreqBins)) {
				T2_Linear += Pxx2[k][i];
			}
			Pxx2_dB[k][i] = 10 * log10(Pxx2[k][i]);
			//									printf(" %d %d %lf\n",k,i,Pxx2_dB[k][i]);
		}
	}

	T2_Linear = T2_Linear / SizeForMean;

	*T2_dB = 10 * log10(T2_Linear) + SysParams->noiseThresh;//average of all the spectrogram in dB + noiseThresh
	return 0;
}

int AvgFilter2(Edge2_Struct *Edge2, int SpectrogramTimeBins, int AvgValue) {
	float filtered_result[SpectrogramTimeBins];
	int SignalLength=SpectrogramTimeBins+AvgValue-1;
	float SignalForAvg[SignalLength];

	int i, m;
	///CHECK AGAIN THE CASE OF PLUS AT MOTION STRUCT 2
	//	for(i=0;i<4;i++)
	//		printf("@@   %d %lf\n",i,Edge2->PrevLastFiftyPrecent[i] );

	for (i = 0; i < SpectrogramTimeBins; i++) {
		filtered_result[i] = 0;
		if (i < AvgValue) {
			for (m = 0; m < AvgValue - i - 1; m++) {//sum on last AvgValue-i bins in current
				filtered_result[i] += Edge2->PrevLastFiftyPrecent[m + i];
			}
			for (m = 0; m <= i; m++) {	//sum on first i bins in current
				filtered_result[i] += Edge2->FiftyPrecent[m];
			}
		} else {
			for (m = 0; m < AvgValue; m++) {
				filtered_result[i] += Edge2->FiftyPrecent[i - m];
			}
		}

		Edge2->FiftyPrecent_Filtered[i] = filtered_result[i] / AvgValue;
		//		printf("### %d %lf\n", i, Edge2->FiftyPrecent_Filtered[i]);

	}

//	for (i = 0; i < SignalLength; i++) {//contacte the end of the previous curve
//	if(i<AvgValue)
//		SignalForAvg[i]=Edge2->PrevLastFiftyPrecent[i];
//
//	else
//			SignalForAvg[i]=Edge2->FiftyPrecent[i-AvgValue];
//	}
//
	return 0;
}

int CalcCurves2(float* Pxx2[], float* Pxx2_dB[], float *T2_dB,
		Edge2_Struct *Edge2, SysParams_Struct* SysParams) {	//CalcRedStarsCurve3 in Matlab
	int k, i, MaxPxx2_dB_idx, idx_FreqAboveThresh, idx_FreqBins_Spectogram = 0,
			AllFreqAboveThresh[SysParams->SpectrogramFreqBins],
			Freq_50_PrecentIdx;
	int MedianType;
	int AvgValue = 5;
	float freq_increment = 1.256281407035176;	//=(Fs/LengthOfDFT);
	float SpectrogramFreqBinsInHz[SysParams->DFTLengthForSpectrogram],
	MaxPxx2_dB, Energy_50_Precent;
	for (float i = 0; i < SysParams->SpectrogramFreqBins; i++) {//create the frequency vector inHz
		SpectrogramFreqBinsInHz[idx_FreqBins_Spectogram] = i * freq_increment;
		idx_FreqBins_Spectogram += 1;
	}
	//	memset(Edge2->Peak,0,174*sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions

	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {

		idx_FreqAboveThresh = 0;
		MaxPxx2_dB = Pxx2_dB[SysParams->Fmin_bin][i];
		MaxPxx2_dB_idx = SysParams->Fmin_bin;
		//		MaxPxx2_dB=Fmin_bin;
		//		SumEnergy=0;
		for (k = SysParams->Fmin_bin;
				k < SysParams->SpectrogramFreqBins - SysParams->noiseFreqBins;
				k++) {//do the work in frequency range [Fmin_bin,SpectrogramFreqBins-noiseFreqBins]
			if (Pxx2_dB[k][i] > *T2_dB) {
				AllFreqAboveThresh[idx_FreqAboveThresh] = k;
				//				SumEnergy+=Pxx2[k][i];
				idx_FreqAboveThresh += 1;
			}

			if (Pxx2_dB[k][i] > MaxPxx2_dB) {		//find MaxPxx2_dB
				MaxPxx2_dB = Pxx2_dB[k][i];
				MaxPxx2_dB_idx = k;
			}
		}

		if (idx_FreqAboveThresh > 0) {//there were frequencies above the threshold
			Edge2->maxFreqIdxs[i] = AllFreqAboveThresh[idx_FreqAboveThresh - 1];
			Edge2->maxFreqEnergy[i] =
					Pxx2_dB[AllFreqAboveThresh[idx_FreqAboveThresh - 1]][i];
			Edge2->PeakEnergy[i] = MaxPxx2_dB;
			Edge2->maxPeakIdxs[i] = MaxPxx2_dB_idx;
			Edge2->Fmax[i] =
					SpectrogramFreqBinsInHz[AllFreqAboveThresh[idx_FreqAboveThresh
															   - 1]];
			//			Edge2->SumEnergy[i]=SumEnergy;
			if (SysParams->truncate == 0)//in this case there is zero-pad at the beginning
				Edge2->Peak[i + SysParams->MedianValue / 2] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx];//first MedianValue/2(=10) bins are 0 for the inital sates of the median filter
			else
				//there is no-zero padding
				Edge2->Peak[i] = SpectrogramFreqBinsInHz[MaxPxx2_dB_idx];
		} else {		//take defaultes
			Edge2->maxFreqIdxs[i] = SysParams->Fmin_bin;
			Edge2->maxFreqEnergy[i] = Pxx2_dB[SysParams->Fmin_bin][i];
			Edge2->PeakEnergy[i] = Pxx2_dB[SysParams->Fmin_bin][i];
			Edge2->maxPeakIdxs[i] = SysParams->Fmin_bin;
			Edge2->Fmax[i] = SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			//			Edge2->SumEnergy[i]=0;
			if (SysParams->truncate == 0)//in this case there is zero-pad at the beginning
				Edge2->Peak[i + SysParams->MedianValue / 2] =
						SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			else
				//there is no-zero padding
				Edge2->Peak[i] = SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
		}

		Energy_50_Precent = (Edge2->PeakEnergy[i] + Edge2->maxFreqEnergy[i])
																																												/ 2;

		for (k = Edge2->maxFreqIdxs[i]; k > Edge2->maxPeakIdxs[i]; k--) {//search backwards the first freq index that passes Energy_50_Precent
			if (Pxx2_dB[k][i] > Energy_50_Precent) {
				break;
			}
		}
		Freq_50_PrecentIdx = k;
		//****think about the case: if ~isempty(Freq_50_PrecentIdx)*****
		Edge2->FiftyPrecent[i] = SpectrogramFreqBinsInHz[Freq_50_PrecentIdx];
		Edge2->FiftyPrecentIdxs[i] = Freq_50_PrecentIdx;

		//							 printf("!!) %d %lf\n",i,Edge2->FiftyPrecent[i] );

	}
	//	MedianFilter2(Edge2, SysParams->SpectrogramTimeBins, SysParams->MedianValue,
	//			SysParams->truncate);//median filter on the peak(black) curve with MedianValue=20
	MedianType=0;//ZERPOAD
	MedianFilter(Edge2, SysParams,
			MedianType);//median filter on the peak(black) curve with MedianValue=20, MedianType=0 ZERPOAD,MedianType=1 TRUNCATE
	for (i = 0; i < SysParams->SpectrogramTimeBins; i++) {

		if (fabs(Edge2->Peak_Filtered[i])
				<= fabs(SpectrogramFreqBinsInHz[SysParams->Fmin_bin])) {// if the peak has no energy, give the 50% also the same energy
			Edge2->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[SysParams->Fmin_bin];
			Edge2->FiftyPrecentIdxs[i] = SysParams->Fmin_bin;
		}

	}
	AvgFilter2(Edge2, SysParams->SpectrogramTimeBins, AvgValue);//Average filter of the FiftyPrecent(blue) curve with AvgValue=5

	return 0;
}

int MaxOfArr(float *Arr, int *MaxIdx, float *MaxValue, int Length) {
	int i;
	*MaxValue = 0;
	for (i = 0; i < Length; i++) {
		if (Arr[i] > *MaxValue) {
			*MaxValue = Arr[i];
			*MaxIdx = i;
		}
	}
	return 0;
}

void MemoryAllocation(Edge2_Struct* Edge2_1, Edge2_Struct* Edge2_2,
		Edge2_Struct* Edge2_Plus_1, Edge2_Struct* Edge2_Minus_1,
		Edge2_Struct* Edge2_Plus_2, Edge2_Struct* Edge2_Minus_2,
		SysParams_Struct *SysParams) {

	Edge2_1->FiftyPrecent = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_1->FiftyPrecent_Filtered = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_1->Fmax = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_1->PeakEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_1->PeakEnergy_PM = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_1->maxFreqEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_1->FreqBins = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_1->maxFreqIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_1->maxPeakIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_1->FiftyPrecentIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_1->PrevLastFiftyPrecent = (float *) malloc(
			(SysParams->AvgValue - 1) * sizeof(float));

	if (SysParams->truncate == 1) {	// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
		Edge2_1->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		Edge2_1->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_1->Peak,0,(SysParams->SpectrogramTimeBins+SysParams->MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_1->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		Edge2_1->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		memset(Edge2_1->Peak, 0,
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Plus_1->FiftyPrecent = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->FiftyPrecent_Filtered = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->Fmax = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->PeakEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->PeakEnergy_PM = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->SumEnergy_Post = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->maxFreqEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->FreqBins = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->maxFreqIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->maxPeakIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->FiftyPrecentIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->PrevLastFiftyPrecent = (float *) malloc(
			(SysParams->AvgValue - 1) * sizeof(float));
	Edge2_Plus_1->T1_t = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Plus_1->SumEnergy_Post, 0,
			(SysParams->SpectrogramTimeBins) * sizeof(float));

	if (SysParams->truncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_1->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		Edge2_Plus_1->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Plus_1->Peak,0,(SysParams->SpectrogramTimeBins+SysParams->MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_1->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		Edge2_Plus_1->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Plus_1->Peak, 0,
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Minus_1->FiftyPrecent = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->FiftyPrecent_Filtered = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->Fmax = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->PeakEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->PeakEnergy_PM = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->SumEnergy_Post = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->maxFreqEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->FreqBins = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->maxFreqIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->maxPeakIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->FiftyPrecentIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->PrevLastFiftyPrecent = (float *) malloc(
			(SysParams->AvgValue - 1) * sizeof(float));
	Edge2_Minus_1->T1_t = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Minus_1->SumEnergy_Post, 0,
			(SysParams->SpectrogramTimeBins) * sizeof(float));

	if (SysParams->truncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_1->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		Edge2_Minus_1->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Minus_1->Peak,0,(SysParams->SpectrogramTimeBins+SysParams->MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the initial conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_1->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		Edge2_Minus_1->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Minus_1->Peak, 0,
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_2->FiftyPrecent = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_2->FiftyPrecent_Filtered = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_2->Fmax = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_2->PeakEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_2->PeakEnergy_PM = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_2->maxFreqEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_2->FreqBins = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_2->maxFreqIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_2->maxPeakIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_2->FiftyPrecentIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_2->PrevLastFiftyPrecent = (float *) malloc(
			(SysParams->AvgValue - 1) * sizeof(float));

	if (SysParams->truncate == 1) {	// Computes medians of smaller segments as it reaches the signal edges-> no zeropadding at the end and begining
		Edge2_2->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		Edge2_2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_1->Peak,0,(SysParams->SpectrogramTimeBins+SysParams->MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints->
		Edge2_2->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		Edge2_2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		memset(Edge2_2->Peak, 0,
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Plus_2->FiftyPrecent = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->FiftyPrecent_Filtered = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->Fmax = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->PeakEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->PeakEnergy_PM = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->SumEnergy_Post = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->maxFreqEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->FreqBins = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->maxFreqIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->maxPeakIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->FiftyPrecentIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->PrevLastFiftyPrecent = (float *) malloc(
			(SysParams->AvgValue - 1) * sizeof(float));
	Edge2_Plus_2->T1_t = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Plus_2->SumEnergy_Post, 0,
			(SysParams->SpectrogramTimeBins) * sizeof(float));

	if (SysParams->truncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_2->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		Edge2_Plus_2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_2->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		Edge2_Plus_2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Plus_2->Peak, 0,
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Minus_2->FiftyPrecent = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->FiftyPrecent_Filtered = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->Fmax = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->PeakEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->PeakEnergy_PM = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->SumEnergy_Post = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->maxFreqEnergy = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->FreqBins = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->maxFreqIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->maxPeakIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->FiftyPrecentIdxs = (int *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->PrevLastFiftyPrecent = (float *) malloc(
			(SysParams->AvgValue - 1) * sizeof(float));
	Edge2_Minus_2->T1_t = (float *) malloc(
			SysParams->SpectrogramTimeBins * sizeof(float));
	memset(Edge2_Minus_2->SumEnergy_Post, 0,
			(SysParams->SpectrogramTimeBins) * sizeof(float));

	if (SysParams->truncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_2->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		Edge2_Minus_2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Minus_1.Peak,0,(SysParams->SpectrogramTimeBins+SysParams->MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the initial conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_2->Peak = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		Edge2_Minus_2->Peak_Filtered = (float *) malloc(
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Minus_2->Peak, 0,
				(SysParams->SpectrogramTimeBins + SysParams->MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

}

void FreeMemory(Edge2_Struct* Edge2_1, Edge2_Struct* Edge2_2,
		Edge2_Struct* Edge2_Plus_1, Edge2_Struct* Edge2_Minus_1,
		Edge2_Struct* Edge2_Plus_2, Edge2_Struct* Edge2_Minus_2,
		Motion_Struct2 *UnitedMotionStruct) {

	free(Edge2_1->FiftyPrecent);
	free(Edge2_1->FiftyPrecent_Filtered);
	free(Edge2_1->Fmax);
	free(Edge2_1->PeakEnergy);
	free(Edge2_1->PeakEnergy_PM);
	free(Edge2_1->maxFreqEnergy);
	free(Edge2_1->FreqBins);
	free(Edge2_1->maxPeakIdxs);
	free(Edge2_1->maxFreqIdxs);
	free(Edge2_1->FiftyPrecentIdxs);
	free(Edge2_1->PrevLastFiftyPrecent);
	free(Edge2_1->Peak);
	free(Edge2_1->Peak_Filtered);

	free(Edge2_Plus_1->FiftyPrecent);
	free(Edge2_Plus_1->FiftyPrecent_Filtered);
	free(Edge2_Plus_1->Fmax);
	free(Edge2_Plus_1->PeakEnergy);
	free(Edge2_Plus_1->PeakEnergy_PM);
	free(Edge2_Plus_1->maxFreqEnergy);
	free(Edge2_Plus_1->FreqBins);
	free(Edge2_Plus_1->maxPeakIdxs);
	free(Edge2_Plus_1->maxFreqIdxs);
	free(Edge2_Plus_1->FiftyPrecentIdxs);
	free(Edge2_Plus_1->PrevLastFiftyPrecent);
	free(Edge2_Plus_1->Peak);
	free(Edge2_Plus_1->Peak_Filtered);

	free(Edge2_Minus_1->FiftyPrecent);
	free(Edge2_Minus_1->FiftyPrecent_Filtered);
	free(Edge2_Minus_1->Fmax);
	free(Edge2_Minus_1->PeakEnergy);
	free(Edge2_Minus_1->PeakEnergy_PM);
	free(Edge2_Minus_1->maxFreqEnergy);
	free(Edge2_Minus_1->FreqBins);
	free(Edge2_Minus_1->maxPeakIdxs);
	free(Edge2_Minus_1->maxFreqIdxs);
	free(Edge2_Minus_1->FiftyPrecentIdxs);
	free(Edge2_Minus_1->PrevLastFiftyPrecent);
	free(Edge2_Minus_1->Peak);
	free(Edge2_Minus_1->Peak_Filtered);

	free(Edge2_2->FiftyPrecent);
	free(Edge2_2->FiftyPrecent_Filtered);
	free(Edge2_2->Fmax);
	free(Edge2_2->PeakEnergy);
	free(Edge2_2->PeakEnergy_PM);
	free(Edge2_2->maxFreqEnergy);
	free(Edge2_2->FreqBins);
	free(Edge2_2->maxPeakIdxs);
	free(Edge2_2->maxFreqIdxs);
	free(Edge2_2->FiftyPrecentIdxs);
	free(Edge2_2->PrevLastFiftyPrecent);
	free(Edge2_2->Peak);
	free(Edge2_2->Peak_Filtered);

	free(Edge2_Plus_2->FiftyPrecent);
	free(Edge2_Plus_2->FiftyPrecent_Filtered);
	free(Edge2_Plus_2->Fmax);
	free(Edge2_Plus_2->PeakEnergy);
	free(Edge2_Plus_2->PeakEnergy_PM);
	free(Edge2_Plus_2->maxFreqEnergy);
	free(Edge2_Plus_2->FreqBins);
	free(Edge2_Plus_2->maxPeakIdxs);
	free(Edge2_Plus_2->maxFreqIdxs);
	free(Edge2_Plus_2->FiftyPrecentIdxs);
	free(Edge2_Plus_2->PrevLastFiftyPrecent);
	free(Edge2_Plus_2->Peak);
	free(Edge2_Plus_2->Peak_Filtered);

	free(Edge2_Minus_2->FiftyPrecent);
	free(Edge2_Minus_2->FiftyPrecent_Filtered);
	free(Edge2_Minus_2->Fmax);
	free(Edge2_Minus_2->PeakEnergy);
	free(Edge2_Minus_2->PeakEnergy_PM);
	free(Edge2_Minus_2->maxFreqEnergy);
	free(Edge2_Minus_2->FreqBins);
	free(Edge2_Minus_2->maxPeakIdxs);
	free(Edge2_Minus_2->maxFreqIdxs);
	free(Edge2_Minus_2->FiftyPrecentIdxs);
	free(Edge2_Minus_2->PrevLastFiftyPrecent);
	free(Edge2_Minus_2->Peak);
	free(Edge2_Minus_2->Peak_Filtered);

	free(UnitedMotionStruct->Edge2->Fmax);
	free(UnitedMotionStruct->Edge2_Plus->Fmax);
	free(UnitedMotionStruct->Edge2_Minus->Fmax);

	free(UnitedMotionStruct->Edge2->Peak_Filtered);
	free(UnitedMotionStruct->Edge2_Plus->Peak_Filtered);
	free(UnitedMotionStruct->Edge2_Minus->Peak_Filtered);

	free(UnitedMotionStruct->Edge2->FiftyPrecent_Filtered);
	free(UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered);
	free(UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered);

	free(UnitedMotionStruct->Edge2_Plus->SumEnergy_Post);
	free(UnitedMotionStruct->Edge2_Minus->SumEnergy_Post);

	free(UnitedMotionStruct->Edge2_Plus->T1_t);
	free(UnitedMotionStruct->Edge2_Minus->T1_t);
}

//more features
//		if(MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]>minfreq){//calculate mean for fifty_precent_filtered without Fmin bins
//			Edge2_noDC_Favg_Plus+=MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];
//			num_of_Edge2_noDC_Favg_Plus+=1;
//		}

//if (num_of_Edge2_noDC_Favg_Plus==0)//equivalent to if isempty(Edge2_noDC) in matlab
//	num_of_Edge2_noDC_Favg_Plus=minfreq;
//Edge2_noDC_Favg_Plus=Edge2_noDC_Favg_Plus/(num_of_Edge2_noDC_Favg_Plus);

//	if(Edge2_noDC_Favg_Plus!=0)//otherwise Edge2_PeakToAvg_Plus stay 0
//		Features_Plus.Edge2_50Precent_PeakToAvg=Edge2_AvgTopFive_Plus/Edge2_noDC_Favg_Plus;
