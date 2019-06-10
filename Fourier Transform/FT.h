#pragma once
#include <iostream>
#include <math.h>
#include <algorithm>
#include <complex>
#include <vector>
class FT
{
public:
	FT();
	~FT();
	void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void FastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void FFT(std::vector<std::complex<double>>&, int h, int u, int v);

	void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseFFT(std::vector<std::complex<double>> &x, int h);

	void LowpassFilter(int** inputImange, int** OutputImage, double ** pFreqReal, double ** pFreqImag, int h, int w, int Cutoff, int n);
	void HighpassFilter(int** inputImange, int** OutputImage, double ** pFreqReal, double ** pFreqImag, int h, int w, int Cutoff, int n);

private:

};



