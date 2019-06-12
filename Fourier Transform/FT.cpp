#include "FT.h"

FT::FT()
{
}

FT::~FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
	}
	for (int forzero_i = 0; forzero_i<M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j<N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // �ť߸��W�v�}�C
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
			//�s�U�ϳť߸���ƻP��Ƴ���
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	// h = w �����O�西
	int N = h;

	// initialize FreqImage
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			FreqReal[i][j] = InputImage[i][j];
			FreqImag[i][j] = 0;
		}
	}

	std::vector<std::complex<double>> Temp;
	Temp.clear();

	// Row FFT
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			std::complex<double> temp(FreqReal[i][j], FreqImag[i][j]);
			Temp.push_back(temp);
		}
		FFT(Temp, h, 0, 0);
		for (int j = 0; j < w; ++j)
		{
			FreqReal[i][j] = Temp[j].real();
			FreqImag[i][j] = Temp[j].imag();
		}
		Temp.clear();
	}

	// Column FFT
	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < h; ++j)
		{
			std::complex<double> temp(FreqReal[j][i], FreqImag[j][i]);
			Temp.push_back(temp);
		}
		FFT(Temp, h, 0, 0);
		for (int j = 0; j < w; ++j)
		{
			FreqReal[j][i] = Temp[j].real();
			FreqImag[j][i] = Temp[j].imag();
		}
		Temp.clear();
	}

	// Output Image
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			OutputImage[i][j] = sqrt(pow(FreqReal[j][i] / N, 2.0f) + pow(FreqImag[j][i] / N, 2.0f));
		}
	}
}	

void FT::FFT(std::vector<std::complex<double>> &x, int h, int u, int v)
{
	// h = w �����O�西
	int N = h;

	// Bit-reversal
	for (int i = 1, j = 0; i < N; ++i)
	{
		for (int k = h >> 1; !((j ^= k)&k); k >>= 1);
		if (i > j)
			swap(x[i], x[j]);
	}

	// FFT
	for (int k = 2; k <= N; k *= 2)
	{
		double omega = -2.0f * 3.14159f / (double)k;
		std::complex<double> dtheta(cos(omega), sin(omega));
		// �Ck�Ӱ��@��FFT
		for (int j = 0; j < N; j += k)
		{
			// �ek/2�ӻP��k/2���T����ƭȫ�n���
			// �]������٪��@�_��
			std::complex<double> theta(1, 0);
			for (int i = j; i < j + k / 2; i++)
			{
				std::complex<double>yeven = x[i];
				std::complex<double>yodd = x[i + k / 2] * theta;
				x[i] = yeven + yodd;
				x[i + k / 2] = yeven - yodd;
				theta *= dtheta;
			}
		}
	}
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	// h = w �����O�西
	int N = h;
	std::vector<std::complex<double>> row;
	row.clear();

	// Row IFFT
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			std::complex<double> temp(FreqReal[i][j], FreqImag[i][j]);
			row.push_back(temp);
		}
		InverseFFT(row, h);
		for (int j = 0; j < w; ++j)
		{
			FreqReal[i][j] = row[j].real() / N;
			FreqImag[i][j] = row[j].imag() / N;
		}
		row.clear();
	}

	// Column IFFT
	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < h; ++j)
		{
			std::complex<double> temp(FreqReal[j][i], FreqImag[j][i]);
			row.push_back(temp);
		}
		InverseFFT(row, h);
		for (int j = 0; j < h; ++j)
		{
			FreqReal[j][i] = row[j].real() / N;
			FreqImag[j][i] = row[j].imag() / N;
		}
		row.clear();
	}

	// Output Image
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			OutputImage[i][j] = sqrt(pow(FreqReal[i][j], 2.0f) + pow(FreqImag[i][j], 2.0f));
		}
	}
}

void FT::InverseFFT(std::vector<std::complex<double>> &x, int h)
{
	// h = w �����O�西
	int N = h;

	// Bit-reversal
	for (int i = 1, j = 0; i < N; ++i)
	{
		for (int k = h >> 1; !((j ^= k)&k); k >>= 1);
		if (i > j)
			swap(x[i], x[j]);
	}

	// IFFT
	for (int k = 2; k <= N; k *= 2)
	{
		double omega = 2.0f * 3.14159f / (double)k;
		std::complex<double> dtheta(cos(omega), sin(omega));
		// �Ck�Ӱ��@��IFFT
		for (int j = 0; j < N; j += k)
		{
			// �ek/2�ӻP��k/2���T����ƭȫ�n���
			// �]������٪��@�_��
			std::complex<double> theta(1, 0);
			for (int i = j; i < j + k / 2; i++)
			{
				std::complex<double>yeven = x[i];
				std::complex<double>yodd = x[i + k / 2] * theta;
				x[i] = yeven + yodd;
				x[i + k / 2] = yeven - yodd;
				theta *= dtheta;
			}
		}
	}
}

void FT::LowpassFilter(int** inputImange, int** OutputImage, double ** pFreqReal, double ** pFreqImag, int h, int w, int Cutoff, int n)
{
	// initial filter
	double** filter = new double*[h];
	for (unsigned int i = 0; i < h; ++i)
		filter[i] = new double[w];
	for (unsigned int i = 0; i < h; ++i)
		for (unsigned int j = 0; j < w; ++j)
			filter[i][j] = 0;

	int u, v, Xmiddle, Ymiddle;
	Xmiddle = w / 2;
	Ymiddle = h / 2;
	// Get filter, find new real, image
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			u = j - Xmiddle;
			v = i - Ymiddle;
			filter[i][j] = 1 / (1 + pow(sqrt(u * u + v * v) / Cutoff, 2 * n));
			std::complex<double> FreTemp(pFreqReal[i][j], pFreqImag[i][j]);
			FreTemp *= filter[i][j];
			pFreqReal[i][j] = FreTemp.real();
			pFreqImag[i][j] = FreTemp.imag();
		}
	}
	for (unsigned int i = 0; i < h; ++i)
	{
		for (unsigned int j = 0; j < w; ++j)
		{
			OutputImage[i][j] = sqrt(pow(pFreqReal[j][i], double(2)) + pow(pFreqImag[j][i], double(2))) / h;
		}
	}
	// Delete filter
	delete[] filter;
}

void FT::HighpassFilter(int** inputImange, int** OutputImage, double ** pFreqReal, double ** pFreqImag, int h, int w, int Cutoff, int n)
{
	// initial filter
	double** filter = new double*[h];
	for (unsigned int i = 0; i < h; ++i)
		filter[i] = new double[w];
	for (unsigned int i = 0; i < h; ++i)
		for (unsigned int j = 0; j < w; ++j)
			filter[i][j] = 0;
	int u, v, Xmiddle, Ymiddle;
	Xmiddle = w / 2;
	Ymiddle = h / 2;
	// Get filter, find new real, image
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			u = j - Xmiddle;
			v = i - Ymiddle;
			filter[i][j] = 1 - (1 / (1 + pow(sqrt(u * u + v * v) / Cutoff, 2 * n)));
			std::complex<double> FreTemp(pFreqReal[i][j], pFreqImag[i][j]);
			FreTemp *= filter[i][j];
			pFreqReal[i][j] = FreTemp.real();
			pFreqImag[i][j] = FreTemp.imag();
		}
	}
	for (unsigned int i = 0; i < h; ++i)
	{
		for (unsigned int j = 0; j < w; ++j)
		{
			OutputImage[i][j] = sqrt(pow(pFreqReal[j][i], double(2)) + pow(pFreqImag[j][i], double(2))) / h;
		}
	}
	// Delete filter
	delete[] filter;
}
