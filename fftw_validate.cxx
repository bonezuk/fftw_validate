#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>

#include "fftw3.h"

//-------------------------------------------------------------------------------------------

const double c_PI_D = 3.141592653589793238464264338832795;

//-------------------------------------------------------------------------------------------
// PRNG is Mersenne Twister implementation MT19937
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//-------------------------------------------------------------------------------------------

#define RANDOM_MT19937_LENGTH_N 624
#define RANDOM_MT19937_LENGTH_M 397
#define RANDOM_MT19937_MATRIX_A 0x9908B0DF
#define RANDOM_MT19937_UPPER_MASK 0x80000000
#define RANDOM_MT19937_LOWER_MASK 0x7fffffff

//-------------------------------------------------------------------------------------------

class Random
{
	public:
		Random();
		~Random();

		// generates a random number on [0,0xffffffff]-interval
		uint32_t randomUInt32();
		// generates a random number on [0,0x7fffffff]-interval
		int32_t randomInt31();
		// generates a random number on [0,0xffffffffffffffff]-interval
		uint64_t randomUInt64();
		
		// generates a random number on [0,1]-real-interval
		double randomReal1();
		// generates a random number on [0,1)-real-interval
		double randomReal2();
		// generates a random number on (0,1)-real-interval
		double randomReal3();
		// generates a random number on [0,1) with 53-bit resolution
		double randomReal53();

		void reset();
		
		void seed(int32_t seedValue);
		void seed(uint32_t seedValue);
		void seed(const uint32_t *key,int32_t keyLength);

	private:

		uint32_t m_mt[RANDOM_MT19937_LENGTH_N];
		int32_t m_mti;
		
		void initializePRNG(uint32_t s);
		void initializePRNGArray(const uint32_t *key,int32_t keyLength);
		
		uint32_t generateRandomUInt32();
};

Random::Random()
{
	reset();
}

//-------------------------------------------------------------------------------------------

Random::~Random()
{}

//-------------------------------------------------------------------------------------------

void Random::reset()
{
	m_mti = RANDOM_MT19937_LENGTH_N + 1;
	initializePRNG(5489UL);
}

//-------------------------------------------------------------------------------------------

void Random::seed(int seedValue)
{
	initializePRNG(static_cast<uint32_t>(seedValue));
}

//-------------------------------------------------------------------------------------------

void Random::seed(uint32_t seedValue)
{
	initializePRNG(seedValue);
}

//-------------------------------------------------------------------------------------------

void Random::seed(const uint32_t *key,int keyLength)
{
	initializePRNGArray(key,keyLength);
}

//-------------------------------------------------------------------------------------------

void Random::initializePRNG(uint32_t s)
{
	m_mt[0] = s & 0xffffffff;
	for(m_mti=1;m_mti < RANDOM_MT19937_LENGTH_N;m_mti++)
	{
		m_mt[m_mti] = (1812433253 * (m_mt[m_mti-1] ^ (m_mt[m_mti-1] >> 30)) + m_mti);
		m_mt[m_mti] &= 0xffffffff;
	}
}

//-------------------------------------------------------------------------------------------

void Random::initializePRNGArray(const uint32_t *key,int keyLength)
{
	int i = 1;
	int j = 0;
	int k = (RANDOM_MT19937_LENGTH_N > keyLength) ? RANDOM_MT19937_LENGTH_N : keyLength;
	
	initializePRNG(19650218UL);
	
	while(k)
	{
		m_mt[i] = (m_mt[i] ^ ((m_mt[i-1] ^ (m_mt[i-1] >> 30)) * 1664525)) + key[j] + j;
		m_mt[i] &= 0xffffffff;
		i++;
		j++;
		if(i>=RANDOM_MT19937_LENGTH_N)
		{
			m_mt[0] = m_mt[RANDOM_MT19937_LENGTH_N-1];
			i=1; 
		}
        if(j>=keyLength)
        {
        	j=0;
        }
		k--;	
	}
	
	for(k=RANDOM_MT19937_LENGTH_N-1;k;k--)
	{
		m_mt[i] = (m_mt[i] ^ ((m_mt[i-1] ^ (m_mt[i-1] >> 30)) * 1566083941)) - i;
		m_mt[i] &= 0xffffffff;
		i++;
		if(i>=RANDOM_MT19937_LENGTH_N)
		{
			m_mt[0] = m_mt[RANDOM_MT19937_LENGTH_N-1];
			i=1;
		}
	}
	
	m_mt[0] = 0x80000000;
}
//-------------------------------------------------------------------------------------------

uint32_t Random::generateRandomUInt32()
{
	static uint32_t mag01[2] = {0x00000000, RANDOM_MT19937_MATRIX_A};

	uint32_t y;

	if(m_mti>=RANDOM_MT19937_LENGTH_N)
	{
		int kk;

		if(m_mti==RANDOM_MT19937_LENGTH_N+1)
		{
			initializePRNG(5489UL);
		}

		for(kk=0;kk<RANDOM_MT19937_LENGTH_N-RANDOM_MT19937_LENGTH_M;kk++)
		{
			y = (m_mt[kk] & RANDOM_MT19937_UPPER_MASK) | (m_mt[kk+1] & RANDOM_MT19937_LOWER_MASK);
			m_mt[kk] = m_mt[kk+RANDOM_MT19937_LENGTH_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		
		for(;kk<RANDOM_MT19937_LENGTH_N-1;kk++)
		{
			y = (m_mt[kk] & RANDOM_MT19937_UPPER_MASK) | (m_mt[kk+1] & RANDOM_MT19937_LOWER_MASK);
			m_mt[kk] = m_mt[kk+(RANDOM_MT19937_LENGTH_M-RANDOM_MT19937_LENGTH_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		
		y = (m_mt[RANDOM_MT19937_LENGTH_N-1] & RANDOM_MT19937_UPPER_MASK) | (m_mt[0] & RANDOM_MT19937_LOWER_MASK);
		m_mt[RANDOM_MT19937_LENGTH_N-1] = m_mt[RANDOM_MT19937_LENGTH_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		m_mti = 0;
	}

	y = m_mt[m_mti++];

	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

//-------------------------------------------------------------------------------------------
// generates a random number on [0,0xffffffff]-interval
//-------------------------------------------------------------------------------------------

uint32_t Random::randomUInt32()
{
	return generateRandomUInt32();
}

//-------------------------------------------------------------------------------------------
// generates a random number on [0,0x7fffffff]-interval
//-------------------------------------------------------------------------------------------

int Random::randomInt31()
{
	return static_cast<int>(generateRandomUInt32() >> 1);
}

//-------------------------------------------------------------------------------------------
// generates a random number on [0,0xffffffffffffffff]-interval
//-------------------------------------------------------------------------------------------

uint64_t Random::randomUInt64()
{
	uint64_t a = static_cast<uint64_t>(generateRandomUInt32());
	uint64_t b = static_cast<uint64_t>(generateRandomUInt32());
	a = ((a << 32) & 0xffffffff00000000ULL) | (b & 0x00000000ffffffffULL);
	return a;
}

//-------------------------------------------------------------------------------------------
// generates a random number on [0,1]-real-interval
//-------------------------------------------------------------------------------------------

double Random::randomReal1()
{
	return static_cast<double>(generateRandomUInt32()) * (1.0/4294967295.0);
}

//-------------------------------------------------------------------------------------------
// generates a random number on [0,1)-real-interval
//-------------------------------------------------------------------------------------------

double Random::randomReal2()
{
	return static_cast<double>(generateRandomUInt32()) * (1.0/4294967296.0);
}

//-------------------------------------------------------------------------------------------
// generates a random number on (0,1)-real-interval
//-------------------------------------------------------------------------------------------

double Random::randomReal3()
{
	return (static_cast<double>(generateRandomUInt32()) + 0.5) * (1.0/4294967296.0);
}

//-------------------------------------------------------------------------------------------
// generates a random number on [0,1) with 53-bit resolution
//-------------------------------------------------------------------------------------------

double Random::randomReal53()
{
	uint32_t a = generateRandomUInt32() >> 5;
	uint32_t b = generateRandomUInt32() >> 6;
	return ((static_cast<double>(a) * 67108864.0) + static_cast<double>(b)) * (1.0/9007199254740992.0);
}


//-------------------------------------------------------------------------------------------
// Complex
//-------------------------------------------------------------------------------------------

class Complex
{
	public:
		Complex();
		Complex(const double& rNum,const double& iNum);
		Complex(const Complex & rhs);
		
		const Complex& operator = (const Complex& rhs);
		
		const Complex& operator += (const Complex& b);
		const Complex& operator -= (const Complex& b);
		const Complex& operator *= (const Complex& b);
		
		double& R();
		const double& R() const;
		double& I();
		const double& I() const;
		
		static Complex W(int n,int m,int N);
		static Complex W(int n,int N);
		
	protected:
	
		double m_R;
		double m_I;
		
		void copy(const Complex& rhs);
};

Complex operator + (const Complex& a,const Complex& b);
Complex operator - (const Complex& a,const Complex& b);
Complex operator * (const Complex& a,const Complex& b);

Complex::Complex() : m_R(0.0),
	m_I(0.0)
{}

//-------------------------------------------------------------------------------------------

Complex::Complex(const double& rNum,const double& iNum) : m_R(rNum),
	m_I(iNum)
{}

//-------------------------------------------------------------------------------------------

Complex::Complex(const Complex & rhs) : m_R(0.0),
	m_I(0.0)
{
	copy(rhs);
}

//-------------------------------------------------------------------------------------------

const Complex& Complex::operator = (const Complex& rhs)
{
	if(this!=&rhs)
	{
		copy(rhs);
	}
	return *this;
}

//-------------------------------------------------------------------------------------------

void Complex::copy(const Complex& rhs)
{
	m_R = rhs.m_R;
	m_I = rhs.m_I;
}

//-------------------------------------------------------------------------------------------

const Complex& Complex::operator += (const Complex& b)
{
	m_R += b.m_R;
	m_I += b.m_I;
	return *this;
}

//-------------------------------------------------------------------------------------------

const Complex& Complex::operator -= (const Complex& b)
{
	m_R -= b.m_R;
	m_I -= b.m_I;
	return *this;
}

//-------------------------------------------------------------------------------------------

const Complex& Complex::operator *= (const Complex& b)
{
	double a = (m_R * b.m_R) - (m_I * b.m_I);
	double c = (m_R * b.m_I) + (m_I * b.m_R);
	m_R = a;
	m_I = c;
	return *this;
}

//-------------------------------------------------------------------------------------------

double& Complex::R()
{
	return m_R;
}

//-------------------------------------------------------------------------------------------

const double& Complex::R() const
{
	return m_R;
}

//-------------------------------------------------------------------------------------------

double& Complex::I()
{
	return m_I;
}

//-------------------------------------------------------------------------------------------

const double& Complex::I() const
{
	return m_I;
}

//-------------------------------------------------------------------------------------------

Complex operator + (const Complex& a,const Complex& b)
{
	Complex c(a);
	c += b;
	return c;
}

//-------------------------------------------------------------------------------------------

Complex operator - (const Complex& a,const Complex& b)
{
	Complex c(a);
	c -= b;
	return c;
}

//-------------------------------------------------------------------------------------------

Complex operator * (const Complex& a,const Complex& b)
{
	Complex c(a);
	c *= b;
	return c;
}

//-------------------------------------------------------------------------------------------

Complex Complex::W(int n,int m,int N)
{
	return W(n * m,N);
}

//-------------------------------------------------------------------------------------------

Complex Complex::W(int n,int N)
{
	double angle = (2.0 * c_PI_D * static_cast<double>(n)) / static_cast<double>(N);
	double R = cos(angle);
	double I = 0.0 - sin(angle);
	return Complex(R,I);
}

//-------------------------------------------------------------------------------------------
// FFTRadix2
//-------------------------------------------------------------------------------------------

template <typename V> class FFTRadix2
{
	public:
		FFTRadix2(int N);
		virtual ~FFTRadix2();
		
		V *DFT(V *x);
        V *IDFT(V *x);
		
	protected:
	
		int32_t m_N;
		int32_t *m_reverseIndex;
		
		V **m_coefficient;
		V **m_stackA;
		V **m_stackB;

		void initialize();
		void done();
		
		int32_t noBits(int32_t N) const;
		int32_t getReverseIndex(int32_t index,int32_t noBits) const;

		V complexMultiplyReal(const V *X,const V *Y) const;
		V complexMultiplyImaginary(const V *X,const V *Y) const;

		V *DFT4(V *x,int index,V *X);
		V *DFT8(V *x,int index,bool sFlag);
		V *DFTRecursive(V *x,int index,int N,bool sFlag);
};

//-------------------------------------------------------------------------------------------

template <typename V> FFTRadix2<V>::FFTRadix2(int N) : m_N(N),
	m_reverseIndex(0),
	m_coefficient(0),
	m_stackA(0),
	m_stackB(0)
{
	initialize();
}

//-------------------------------------------------------------------------------------------

template <typename V> FFTRadix2<V>::~FFTRadix2()
{
	done();
}

//-------------------------------------------------------------------------------------------

template <typename V> void FFTRadix2<V>::initialize()
{
	int32_t nBits = noBits(m_N);
	m_reverseIndex = new int32_t [m_N];
	for(int32_t i=0;i<m_N;i++)
	{
		m_reverseIndex[i] = getReverseIndex(i,nBits);
	}
	
	m_coefficient = reinterpret_cast<V **>(malloc((nBits - 2) * sizeof(V *)));
	m_stackA = reinterpret_cast<V **>(malloc((nBits - 2) * sizeof(V *)));
	m_stackB = reinterpret_cast<V **>(malloc((nBits - 3) * sizeof(V *)));

	double halfSquared = sqrt(0.5);
	V *c2 = new V [4];
	c2[0] =  halfSquared;
	c2[1] = -halfSquared;
	c2[2] = -halfSquared;
	c2[3] = -halfSquared;
	m_coefficient[0] = c2;
	m_stackA[0] = new V [32];
	m_stackB[0] = new V [32];
	
	for(int32_t i=4;i<=nBits;i++)
	{
		int32_t M = 1 << i;
		int32_t len = (M >> 1) - 1;
		V *coefficients = new V [len << 1];
		for(int32_t j=1; j<(M >> 1); j++)
		{
			double angle = (2.0 * c_PI_D * static_cast<double>(j)) / static_cast<double>(M);
			coefficients[((j-1)<<1)+0] = static_cast<V>(cos(angle));
			coefficients[((j-1)<<1)+1] = static_cast<V>(0.0 - sin(angle));
		}
		m_coefficient[i-3] = coefficients;
		if(i<nBits)
		{
			m_stackA[i-3] = new V [M<<1];
			m_stackB[i-3] = new V [M<<1];
		}
	}
}

//-------------------------------------------------------------------------------------------

template <typename V> void FFTRadix2<V>::done()
{
	int32_t nBits = noBits(m_N);

	if(m_reverseIndex!=0)
	{
		delete [] m_reverseIndex;
		m_reverseIndex = 0;
	}
	if(m_coefficient!=0)
	{
		for(int32_t i=3;i<=nBits;i++)
		{
			V *c = m_coefficient[i-3];
			delete [] c;
			if(i<nBits)
			{
				V *sA = m_stackA[i-3];
				delete [] sA;
				V *sB = m_stackB[i-3];
				delete [] sB;
			}
		}
		free(m_coefficient);
		m_coefficient = 0;
		free(m_stackA);
		m_stackA = 0;
		free(m_stackB);
		m_stackB = 0;
	}
}

//-------------------------------------------------------------------------------------------

template <typename V> int32_t FFTRadix2<V>::noBits(int32_t N) const
{
	int32_t count = 0;
	
	while(N > 1)
	{
		N >>= 1;
		count++;
	}
	return count;
}

//-------------------------------------------------------------------------------------------

template <typename V> int32_t FFTRadix2<V>::getReverseIndex(int32_t index,int32_t noBits) const
{
	uint32_t y = static_cast<uint32_t>(index), x = 0;
	
	while(noBits>0)
	{
		x = (x << 1) | (y & 0x00000001);
		y >>= 1;
		noBits--;
	}
	return static_cast<int32_t>(x << 1);
}

//-------------------------------------------------------------------------------------------

template <typename V> V FFTRadix2<V>::complexMultiplyReal(const V *X,const V *Y) const
{
	return (X[0] * Y[0]) - (X[1] * Y[1]);
}

//-------------------------------------------------------------------------------------------

template <typename V> V FFTRadix2<V>::complexMultiplyImaginary(const V *X,const V *Y) const
{
	return (X[0] * Y[1]) + (X[1] * Y[0]);
}

//-------------------------------------------------------------------------------------------

template <typename V> V *FFTRadix2<V>::DFT4(V *x,int index,V *X)
{
	int32_t idxA,idxB,idxC,idxD;
	V F1[4],F2[4];

	idxA = m_reverseIndex[index+0];
	idxB = m_reverseIndex[index+1];
	idxC = m_reverseIndex[index+2];
	idxD = m_reverseIndex[index+3];

	F1[0] = x[idxA] + x[idxB];
	F1[1] = x[idxA+1] + x[idxB+1];
	F1[2] = x[idxA] - x[idxB];
	F1[3] = x[idxA+1] - x[idxB+1];

	F2[0] = x[idxC] + x[idxD];
	F2[1] = x[idxC+1] + x[idxD+1];
	F2[2] = x[idxC] - x[idxD];
	F2[3] = x[idxC+1] - x[idxD+1];

	X[0] = F1[0] + F2[0];
	X[1] = F1[1] + F2[1];
	X[2] = F1[2] + F2[3];
	X[3] = F1[3] - F2[2];
	X[4] = F1[0] - F2[0];
	X[5] = F1[1] - F2[1];
	X[6] = F1[2] - F2[3];
	X[7] = F1[3] + F2[2];
		
	return X;
}

//-------------------------------------------------------------------------------------------

template <typename V> V *FFTRadix2<V>::DFT8(V *x,int index,bool sFlag)
{
	V *coefficients = m_coefficient[0];
	V *X = (sFlag) ? m_stackA[0] : m_stackB[0];

	V t;
	V *F1 = DFT4(x,index,&X[16]);
	V *F2 = DFT4(x,index+4,&X[24]);
	
	t = complexMultiplyReal(&coefficients[0],&F2[2]);
	F2[3] = complexMultiplyImaginary(&coefficients[0],&F2[2]);
	F2[2] = t;

	t = complexMultiplyReal(&coefficients[2],&F2[6]);
	F2[7] = complexMultiplyImaginary(&coefficients[2],&F2[6]);
	F2[6] = t;
	
	X[ 0] = F1[0] + F2[0];
	X[ 1] = F1[1] + F2[1];
	X[ 2] = F1[2] + F2[2];
	X[ 3] = F1[3] + F2[3];
	X[ 4] = F1[4] + F2[5];
	X[ 5] = F1[5] - F2[4];
	X[ 6] = F1[6] + F2[6];
	X[ 7] = F1[7] + F2[7];
	
	X[ 8] = F1[0] - F2[0];
	X[ 9] = F1[1] - F2[1];
	X[10] = F1[2] - F2[2];
	X[11] = F1[3] - F2[3];
	X[12] = F1[4] - F2[5];
	X[13] = F1[5] + F2[4];
	X[14] = F1[6] - F2[6];
	X[15] = F1[7] - F2[7];
	
	return X;
}

//-------------------------------------------------------------------------------------------

template <typename V> V *FFTRadix2<V>::DFTRecursive(V *x,int index,int N,bool sFlag)
{
	int32_t i;
	int32_t halfN = N >> 1;
	int32_t bitIndex = noBits(N) - 3;
	V *coefficients = m_coefficient[bitIndex];
	V *X;
	
	if(N==m_N)
	{
		X = new V [N << 1];
	}
	else
	{
		X = (sFlag) ? m_stackA[bitIndex] : m_stackB[bitIndex];
	}
	
	V *Y = &X[N];
	V *F1 = (halfN != 8) ? DFTRecursive(x,index,halfN,true) : DFT8(x,index,true);
	V *F2 = (halfN != 8) ? DFTRecursive(x,index+halfN,halfN,false) : DFT8(x,index+halfN,false);
	
	for(i=1;i<halfN;i++)
	{
		int32_t idxA = (i-1) << 1;
		int32_t idxB = i << 1;
		double t = complexMultiplyReal(&coefficients[idxA],&F2[idxB]);
		F2[idxB+1] = complexMultiplyImaginary(&coefficients[idxA],&F2[idxB]);
		F2[idxB] = t;
	}
	
	for(i=0;i<N;i++)
	{
		X[i] = F1[i] + F2[i];
		Y[i] = F1[i] - F2[i];
	}
	
	return X;
}

//-------------------------------------------------------------------------------------------

template <typename V> V *FFTRadix2<V>::DFT(V *x)
{
	return DFTRecursive(x,0,m_N,false);
}

//-------------------------------------------------------------------------------------------

template <typename V> V *FFTRadix2<V>::IDFT(V *x)
{
    V *X = DFTRecursive(x,0,m_N,false);
    for(int i = 0; i < m_N << 1; i++)
    {
        X[i] /= static_cast<double>(m_N);
    }
    return X;
}

//-------------------------------------------------------------------------------------------

Complex *DFT_N_Full(Complex *x,int N)
{
	Complex *X = new Complex [N];
	
	for(int k=0;k<N;k++)
	{
		Complex c;
		
		for(int n=0;n<N;n++)
		{
			c += x[n] * Complex::W(k,n,N);
		}
		X[k] = c;
	}
	return X;
}

//-------------------------------------------------------------------------------------------

Complex *IDFT_N_Full(Complex *x,int N)
{
    double invN = 1.0 / static_cast<double>(N);
    Complex *X = DFT_N_Full(x, N);
    for(int k = 0; k < N; k++)
    {
        X[k].R() *= invN;
        X[k].I() *= invN;
    }
	return X;
}

//-------------------------------------------------------------------------------------------

void DFT_R2C_Full(int N, double *inR, double *outC)
{
	for(int i = 0; i < N; i++)
	{
		double sumR = 0.0, sumI = 0.0;
		for(int j = 0; j < N; j++)
		{
			sumR += inR[j] * cos((2.0 * c_PI_D * static_cast<double>(i * j)) / static_cast<double>(N));
			sumI -= inR[j] * sin((2.0 * c_PI_D * static_cast<double>(i * j)) / static_cast<double>(N));
		}
		outC[(i << 1) + 0] = sumR;
		outC[(i << 1) + 1] = sumI;
	}
}

//-------------------------------------------------------------------------------------------

void InverseDFT_C2R_Full(int N, double *inC, double *outR)
{
	for(int i = 0; i < N; i++)
	{
		double sum = 0.0;
		for(int j = 0; j < N; j++)
		{
			sum += inC[(j << 1) + 0] * cos((2.0 * c_PI_D * static_cast<double>(i * j)) / static_cast<double>(N));
			sum -= inC[(j << 1) + 1] * sin((2.0 * c_PI_D * static_cast<double>(i * j)) / static_cast<double>(N));
		}
		outR[i] = sum / static_cast<double>(N);
	}
}

//-------------------------------------------------------------------------------------------

void runDFTTest(int size)
{
	const double c_TOLERANCE = 0.00000001;
	int i;
	Random rand;
	rand.seed(0);

	double* inA = reinterpret_cast<double*>(fftw_malloc(size * sizeof(double)));
	assert(inA != NULL);
	fftw_complex* freqA = reinterpret_cast<fftw_complex*>(fftw_malloc(((size / 2) + 1) * sizeof(fftw_complex)));
	assert(freqA != NULL);
	double* outA = reinterpret_cast<double*>(fftw_malloc(size * sizeof(double)));
	assert(outA != NULL);

	double *inB = new double [size];
	double *freqB = new double [size * 2];
	double *outB = new double [size];

	Complex *inC = new Complex [size];
	double *inD = new double [size * 2];

	for(i = 0; i < size; i++)
	{
		inA[i] = rand.randomReal1();
		inB[i] = inA[i];
		inC[i].R() = inA[i];
		inC[i].I() = 0.0;
		inD[(i << 1) + 0] = inA[i];
		inD[(i << 1) + 1] = 0.0;
	}

	fftw_plan planF = fftw_plan_dft_r2c_1d(size, inA, freqA, FFTW_ESTIMATE);
	assert(planF != NULL);
	fftw_plan planI = fftw_plan_dft_c2r_1d(size, freqA, outA, FFTW_ESTIMATE);
	assert(planI != NULL);
	fftw_execute(planF);
	fftw_execute(planI);
	for(i = 0; i < size; i++)
	{
		outA[i] /= static_cast<double>(size);
	}

	DFT_R2C_Full(size, inB, freqB);
	InverseDFT_C2R_Full(size, freqB, outB);

	Complex *freqC = DFT_N_Full(inC, size);
	Complex *outC = IDFT_N_Full(freqC, size);

	FFTRadix2<double> dft(size);
	double *freqD = dft.DFT(inD);
	double *outD = dft.IDFT(freqD);

	char fileName[512];
	snprintf(fileName, 512, "/Users/bonez/Development/tmp/dft_%d.csv", size);
	FILE *fout = fopen(fileName, "w");
	assert(fout != NULL);
	fprintf(fout, "idx,inA,freqA.R,freqA.I,outA,freqB.R,freqB.I,outB,freqC.R,freqC.I,outC,freqD.R,freqD.I,outD,diffB,diffC,diffD\n");
	for(i = 0; i < size / 2; i++)
 	{
		fprintf(fout, "%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,\n", i,
			inA[i], freqA[i][0], freqA[i][1], outA[i],
			freqB[(i << 1) + 0], freqB[(i << 1) + 1], outB[i],
			freqC[i].R(), freqC[i].I(), outC[i].R(),
			freqD[(i << 1) + 0], freqD[(i << 1) + 1], outD[(i << 1) + 0],
			freqA[i][0] - freqB[(i << 1) + 0], freqA[i][0] - freqC[i].R(), freqA[i][0] - freqD[(i << 1) + 0]);
 	}
	fflush(fout);
	fclose(fout);

	delete [] inD;
	delete [] freqD;
	delete [] outD;

	delete [] inC;
	delete [] freqC;
	delete [] outC;

	delete [] inB;
	delete [] freqB;
	delete [] outB;

	fftw_destroy_plan(planF);
	fftw_destroy_plan(planI);
	fftw_free(outA);
	fftw_free(freqA);
	fftw_free(inA);
}

//-------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	runDFTTest(64);
	runDFTTest(128);
	runDFTTest(256);
	runDFTTest(512);
	return 0;
}

//-------------------------------------------------------------------------------------------
