

#ifndef RANDOM_NUMBERS_CPP_H
#define RANDOM_NUMBERS_CPP_H

class RandomNumbers
{
public:
	RandomNumbers():seed31(1002002){}
	RandomNumbers(int seed){ seed31 = seed==0?1:seed; }

	void init(int seed){ seed31 = seed==0?1:seed; }
	int randInt(int min, int max);
	int randIntDifferentThan(int min, int max, int diff);
	double rand01();
	double randNormal(double mean, double std_dev);
	double randBeta(double alpha, double beta);
	double randGamma(double alpha);
	double randBetaJohnk(double alpha, double beta);
	double randBetaMarsagliaTsang(double alpha, double beta);

private:
	long unsigned int next_rand();
	long unsigned int seed31;
};

#endif
