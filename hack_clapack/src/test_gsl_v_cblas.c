#include "lap2gsl.c"
#include <f2c.h>
#include <clapack.h>

void printMatrix(int len, double * matrix)
{
	for(int i=0;i<len;i++)
	{
		for(int j=0; j<len; j++)
		{
			printf("%f ", matrix[i*len+j]);
		}
		printf("\n");
	}
	return;
}

void printVector(int len, double * vector)
{
	for(int i=0; i<len; i++)
	{
		printf("%f ", vector[i]);
	}
	printf("\n");
	return;
}

void printResults(int len, double * matrix, double * vector)
{
	printf("Matrix result of:\n");
	printMatrix(len,matrix);

	printf("\nand a solution of:\n");
	printVector(len,vector);
	return;

}

int main()
{
//test the dposv_ code
double matrixa[16] = {990,  556,  661, 1103,  556,  333,  356,  608,  661,  356,  583,
        773, 1103,  608,  773, 1278};
double vectorb[4] = {1233,3214,7646,1237};
double matrixc[16] = {990,  556,  661, 1103,  556,  333,  356,  608,  661,  356,  583,
        773, 1103,  608,  773, 1278};	
double vectord[4] = {1233,3214,7646,1237};
double expected[4] = {-118.90029333,  163.78255239,   70.08079298,  -16.72004796};

long one =1;
long info = 0;
long num =4;
char side = 'L';

//clapack version
dposv_(&side,&num,&one,matrixa,&num,vectorb,&num,&info);
//gsl versioin
mdposv_(&side,&num,&one,matrixc,&num,vectord,&num,&info);

//print the results from clapack
printf("clapack\n");
printResults(4,matrixa,vectorb);
//print the results from gls
printf("gls\n");
printResults(4,matrixc,vectord);
printf("Expected result:\n");
printVector(num,expected);
printf("\n");


//test the inversion of a matrix

double matrixe[16] = {990,  556,  661, 1103,  556,  333,  356,  608,  661,  356,  583,
        773, 1103,  608,  773, 1278};
double matrixf[16] = {990,  556,  661, 1103,  556,  333,  356,  608,  661,  356,  583,
        773, 1103,  608,  773, 1278};	
double expectedInverseL[16] = {1.01010101e-03,  -1.68653502e-03,  -1.15386450e-04,
         3.64322476e-07,   0.00000000e+00,   3.00300300e-03,
        -1.83373768e-03,  -3.19520029e-04,   0.00000000e+00,
         0.00000000e+00,   1.71526587e-03,  -1.03748084e-03,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         7.82472613e-04};
double expectedInverseU[16] = {1.01010101e-03,  -2.89382701e-19,  -1.13145981e-19,
         4.38747026e-20,  -1.68653502e-03,   3.00300300e-03,
         9.13204407e-19,  -2.90383362e-19,  -1.15386450e-04,
        -1.83373768e-03,   1.71526587e-03,   3.15871823e-19,
         3.64322476e-07,  -3.19520029e-04,  -1.03748084e-03,
         7.82472613e-04};


printf("\nTest dtrtri\n");
//clapack version
dtrtri_("L","N",&num,matrixe,&num,&info);
//gsl version
mdtrtri_("L","N",&num,matrixf,&num,&info);
printf("\nThe return code of gsl dtrtri is %ld\n",info);

//print the results from clapack
printf("\nThe results from clapack are\n");
printMatrix(num,matrixe);
//print the results from gsl
printf("\nThe results from gsl are\n");
printMatrix(num,matrixf);
//expected results
printf("\nThe expected inversion for Lower triangle\n");
printMatrix(num,expectedInverseL);
return 0;
}
