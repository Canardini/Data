#include<Volatility/GaussianFunctional.h>


GFModel::GFModel(double inSkew, double inConvexity, double inLLSkew, double inRRskew, double inRCMS, double inBeta,
	 double inATMBPVol, double inMaturity, double inFWD) :mRCMS(inRCMS), mSkew(inSkew), mBeta(inBeta), mATMBPVol(inATMBPVol), mConvexity(inConvexity),
	mLLSkew(inLLSkew), mRRSkew(inRRskew), mMaturity(inMaturity), mFwd(inFWD), mFWDGrid(5)

{
	mCutOff = 4.0;
	mStdDev = inATMBPVol*sqrt(inMaturity);
	mPositivityEpsilon = mStdDev*0.05;
	mF0 = mFwd;
	mFPOne = mF0 + mStdDev;
	mFMOne = mF0 - mStdDev;
	mFPFour = mF0 + mCutOff*mStdDev;
	mFMFour = mF0 - mCutOff*mStdDev;
	mFMax = fmax(0.4, mF0 + 10.0*mStdDev);
	mFMin = fmin(0.0, mF0 - 5.0*mStdDev) - 1.0;
	mFWDGrid[0] = mFMin;
	mFWDGrid[1] = mFMFour;
	mFWDGrid[2] = mFMOne;
	mFWDGrid[3] = mFPOne;
	mFWDGrid[4] = mFPFour;
	mCenterVol = inATMBPVol / backbone(mFwd);
	mSQRTT = sqrt(mMaturity);

	int lUniformPoints(100);
	double lUniformStep = 1.0 / lUniformPoints;
	mStates.resize(lUniformPoints + 1);
	mLowerBound = -10;
	mUpperBound = 10;
	double lU;
	lUniformStep = (mUpperBound - mLowerBound) / lUniformPoints;
	mStateStep = lUniformStep;

	for (int i = 0; i < lUniformPoints + 1; ++i){
		lU = i*lUniformStep;
		mStates[i] = mLowerBound + lU;;
	}

	mCenterVol = getAlpha(mATMBPVol);

}


double GFModel::getAlpha(double inATMVol){

	double lSigmaNew = 1e100;
	double lSigma0 = inATMVol;
	double lF, lFp;
	int lCompteur(0);
	while ((abs(lSigmaNew - lSigma0) > 1e-10)&(lCompteur<20)){
		lSigmaNew = lSigma0;
		mCenterVol = lSigmaNew;
		this->buildVolFunction();
		lF = this->getNormalVol(mFwd) - inATMVol;
		mCenterVol = lSigmaNew + 1e-6;
		this->buildVolFunction();
		lFp = this->getNormalVol(mFwd) - inATMVol - lF;
		lSigma0 = lSigmaNew - (lF / lFp)*1e-6;
		lCompteur++;
	}
	return lSigma0;
}


void GFModel::buildVolFunction(){


	mSkewTotal = mSkew + mATMBPVol*diffBackbone(mFwd);
	double lSigma0 = mCenterVol;
	double lSigmaPOne = lSigma0 + (mSkewTotal + mConvexity)*mStdDev;
	double lSigmaMOne = lSigma0 + (mSkewTotal - mConvexity)*(-mStdDev);
	double lSigmaMFour = lSigma0 + (mSkewTotal - mConvexity - mLLSkew)*(-mCutOff*mStdDev);
	double lSigmaPFour = lSigma0 + (mSkewTotal + mConvexity + mRRSkew)*(mCutOff*mStdDev);
	double lSigmaPMax = lSigmaPFour + (10 - mCutOff)*mRCMS*mStdDev;

	//Center polynomial (2nd order)
	getQuadraticParams(mFMOne, mF0, mFPOne, lSigmaMOne, lSigma0, lSigmaPOne, mCenterParam0, mCenterParam1, mCenterParam2);

	//Right side 
	getLinearParams(mFPOne, mFPFour, lSigmaPOne, lSigmaPFour, mRightParam0, mRightParam1);
	//Extreme Right (CMS)
	getLinearParams(mFPFour, mFMax, lSigmaPFour, lSigmaPMax, mRightRightParam0, mRightRightParam1);

	//Left side 
	getLinearParams(mFMFour, mFMOne, lSigmaMFour, lSigmaMOne, mLeftParam0, mLeftParam1);

	/// Left Left side
	getTwoPointsCubicParams(mFMin, mFMFour, 0, lSigmaMFour, 0, mLeftParam0, mLeftLeftParam0, mLeftLeftParam1, mLeftLeftParam2, mLeftLeftParam3);

}

double GFModel::getVolatility(double inX){

	if ((inX >= mFMOne)&(inX <= mFPOne)){
		double lLocalVol = mCenterParam0 + mCenterParam1*inX + mCenterParam2*inX*inX;
		lLocalVol = positivity(lLocalVol);
		return lLocalVol;
	}
	else if (inX < mFMOne){
		if (inX < mFMFour){
			double lLocalVol = mLeftLeftParam0 + mLeftLeftParam1*inX + mLeftLeftParam2*inX*inX + mLeftLeftParam3*inX*inX*inX;
			lLocalVol = positivity(lLocalVol);
			return lLocalVol;
		}
		else{
			double lLocalVol = mLeftParam0 + mLeftParam1*inX;
			lLocalVol = positivity(lLocalVol);
			return lLocalVol;

		}
	}
	else{
		if (inX <= mFPFour){
			double lLocalVol = mRightParam0 + mRightParam1*inX;
			lLocalVol = positivity(lLocalVol);
			return lLocalVol;
		}
		else{
			double lLocalVol = mRightRightParam0 + mRightRightParam1*inX;
			lLocalVol = positivity(lLocalVol);
			return lLocalVol;
		}
	}
}




void GFModel::getLinearParams(double inX1, double inX2, double inB1, double inB2, double & outA0, double & outA1){
	double lRatio1 = (inB1 - inB2) / (inX1 - inX2);
	outA1 = lRatio1;
	outA0 = inB1 - outA1*inX1;
}


void GFModel::getQuadraticParams(double inX1, double inX2, double inX3, double inB1, double inB2, double inB3, double & outA0, double & outA1, double &outA2){
	double lRatio1 = (inB1 - inB2) / (inX1 - inX2);
	double lRatio2 = (inB2 - inB3) / (inX2 - inX3);
	outA2 = (lRatio2 - lRatio1) / (inX3 - inX1);
	outA1 = lRatio2 - outA2*(inX2 + inX3);
	outA0 = inB1 - outA1*inX1 - outA2*inX1*inX1;
}


void GFModel::getTwoPointsCubicParams(double inX1, double inX2, double inB1, double inB2, double inB3,
	double inB4, double & outA0, double & outA1, double &outA2, double &outA3){
	double lDenominator = 0.5*(inX1*inX1 + inX2*inX2) - inX1*inX2;
	double lRatio1 = (inB1 - inB2) / (inX1 - inX2);
	double lRatio2 = (inB3 - inB4) / (inX1 - inX2);

	outA3 = (inB4 - lRatio1 + lRatio2*(0.5*(inX1 - inX2))) / lDenominator;
	outA2 = 0.5*(lRatio2 - 3 * outA3*(inX1 + inX2));
	outA1 = lRatio1 - 0.5*lRatio2*(inX1 + inX2) + outA3*(0.5*(inX1*inX1 + inX2*inX2) + 2 * inX1*inX2);
	outA0 = inB1 - outA1*inX1 - outA2*inX1*inX1 - outA3*inX1*inX1*inX1;
}

double GFModel::backbone(double inX){
	if (inX <= 0){
		double lU = inX / 0.01;
		return pow((1 - lU + lU* lU), -mBeta);
	}
	else if (inX <= 0.08){
		double lU = inX / 0.01;
		double lU8 = inX / 0.08;
		return pow((1 + lU - 10 * pow(lU8, 3) + 15 * pow(lU8, 4) - 6 * pow(lU8, 5)), -mBeta);
	}
	else{
		double lU = inX / 0.01;
		return pow(lU, mBeta);
	}
}


double GFModel::diffBackbone(double inX){
	///return (B(F')/B(F)-1)/epsilon
	double lEpsilon(1e-6);
	double lRes = backbone(inX + lEpsilon) / backbone(inX) - 1.0;

	return lRes / lEpsilon;
}

double GFModel::positivity(double inX){
	if (inX >= mPositivityEpsilon){
		return inX;
	}
	else{
		return (inX - 2 * mPositivityEpsilon)*mPositivityEpsilon / (2 * inX - 3 * mPositivityEpsilon);
	}
}
 
 
 
double GFModel::nextR4(double inStep, double inFactor, double inF){
	double k1 = inFactor*this->getVolatility(inF);
	double k2 = inFactor*this->getVolatility(inF + inStep*0.5*k1);
	double k3 = inFactor*this->getVolatility(inF + inStep*0.5*k2);
	double k4 = inFactor*this->getVolatility(inF + inStep*k3);

	double lRes = inF + inStep*(k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
	return lRes;
}

 

void GFModel::getODEFWD2(double inFactor, std::vector<double> & loutFWDS){

	int lSize(loutFWDS.size() / 2);
	loutFWDS[lSize] = mFwd;

	for (int i = lSize; i < loutFWDS.size() - 1; ++i){
		loutFWDS[i + 1] = nextR4(mStateStep, inFactor, loutFWDS[i]);
	}
	for (int i = lSize; i >0; --i){
		loutFWDS[i - 1] = nextR4(-mStateStep, inFactor, loutFWDS[i]);
	}
}


void GFModel::getLeftQuadraticTwoPointsParams(double inX1, double inX2, double inB1, double inB2, double inB3, double & outA, double & outB, double &outC){
	double lRatio = (inB1 - inB2) / (inX1 - inX2);
	outA = -lRatio / (inX1 - inX2) + inB3 / (inX1 - inX2);
	outB = inB3 - 2 * outA*inX1;
	outC = inB1 - outA*inX1*inX1 - outB*inX1;
}

void GFModel::getRightQuadraticTwoPointsParams(double inX1, double inX2, double inB1, double inB2, double inB3, double & outA, double & outB, double &outC){
	double lRatio = (inB1 - inB2) / (inX1 - inX2);
	outA = lRatio / (inX1 - inX2) - inB3 / (inX1 - inX2);
	outB = inB3 - 2 * outA*inX2;
	outC = inB1 - outA*inX1*inX1 - outB*inX1;
}

double GFModel::getPrice(double inStrike){
	
	std::vector<double> lFWDs(mStates.size());
	int lRefIndex(mStates.size() / 2);
	bool lFound(false);
	double lRefIncreas = -1e100;
	bool isIncreased(true);

	getODEFWD2(mSQRTT, lFWDs);
	for (int i = 0; i <mStates.size(); ++i){
		if ((lFWDs[i]>inStrike)&(!lFound)){
			lFound = true;
			lRefIndex = i;
		}
		if ((lFWDs[i] < lRefIncreas) &isIncreased){
			isIncreased = false;
		}
		lRefIncreas = lFWDs[i];
	}
	if (!isIncreased){
		std::sort(lFWDs.begin(), lFWDs.end());
	}

	MonotonicCubicSpline2 oSpline(mStates, lFWDs);

	double lResIntegral = oSpline.getTotalIntegral();
	double lEpsilon = mFwd - lResIntegral;
	double lEffectiveStrike = inStrike - lEpsilon;
	int lLocateIndex;
	double lStartingState;

	if (lFWDs[lFWDs.size() - 1] <= lEffectiveStrike){
		lLocateIndex = lFWDs.size() - 1;
	}
	else if (lFWDs[0] >= lEffectiveStrike){
		lLocateIndex = -1;
	}
	else{

		lLocateIndex = oSpline.locateY(lEffectiveStrike);
		int lCompteur(1);
		double a, b, c;
		a = mStates[lLocateIndex];
		b = mStates[lLocateIndex + 1];
		double lfa, lfb, lfc;
		lfa = oSpline.getY(a) - lEffectiveStrike;
		lfb = oSpline.getY(b) - lEffectiveStrike;
		while (lCompteur <= 50){
			c = 0.5*(a + b);
			lfc = oSpline.getY(c) - lEffectiveStrike;
			if ((lfc == 0) | (0.5*(b - a) < 1e-10)){
				break;
			}
			lCompteur++;
			if (lfa*lfc > 0){
				a = c;
				lfa = lfc;
			}
			else{
				b = c;
				lfb = lfc;
			}
		}
		lStartingState = fmin(fmax(c, mLowerBound), mUpperBound);
		lLocateIndex = oSpline.locate(lStartingState);
	}

	double lCallPrice;
	if (lLocateIndex == -1){
		lCallPrice = mFwd - inStrike;
	}
	else if (lLocateIndex == mStates.size() - 1){
		lCallPrice = 0.0;
	}
	else{
		double lIntegral = oSpline.getNormalIntegral(lStartingState, lLocateIndex);
		for (int i = lLocateIndex + 1; i < mStates.size() - 1; ++i){
			lIntegral += oSpline.getIntegrals(i);
		}
		lCallPrice = lIntegral - lEffectiveStrike*Normalcdf(-lStartingState);
	}

	return lCallPrice;
}


double GFModel::getNormalVol(double  inStrike){
	double lBpsVol;
	double lCallPrice = this->getPrice(inStrike);
	lBpsVol = PriceToBpsVol(lCallPrice, mMaturity, mFwd, inStrike, "C");
	return lBpsVol;
}