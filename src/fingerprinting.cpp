#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppGSL)]]
using namespace Rcpp;

gsl_rng *rng;

typedef struct {
	double gof1;
	double gof2;
	double *w;
} type_try;

// Comparison function to sort with
int compare(const void *p1, const void *p2)
{
	type_try *elem1 = (type_try*)p1;    
	type_try *elem2 = (type_try*)p2;

	if ( elem1->gof2 < elem2->gof2 )
	  return 1;
	else if ( elem1->gof2 > elem2->gof2 )
	  return -1;
	else
	  return 0;
}

// Comparison function to sort with
int compare2(const void *p1, const void *p2)
{
	double *n1 = (double*)p1;    
	double *n2 = (double*)p2;

	if ( *n1 > *n2 )
	  return 1;
	else if ( *n1 < *n2 )
	  return -1;
	else
	  return 0;
}

double correct(double avg, double dev, int n)
{
	if( n > 1 )
	{
		// return avg + gsl_ran_tdist(rng, (double)(n-1)) * dev;
		return avg + gsl_ran_tdist(rng, (double)(n-1)) * dev / sqrt((double)n);
	}
	return avg;
}

//' unmix_c
//'
//' @param sources Data frame containing sediment source samples
//' @param samples Data frame containing mixture samples
//' @param trials Number of samples in each hypercube dimension
//' @param iter Iterations in the source variability analysis
//' @param seed Seed for the random number generator
//' @return Data frame containing the relative contribution of the sediment sources for each sediment mixture and iterations
// [[Rcpp::export]]
Rcpp::DataFrame unmix_c(SEXP sources, SEXP samples, int trials=100, int iter=100, int seed=69512)
{
	// Construct the data.frame object
	Rcpp::DataFrame dfs = Rcpp::DataFrame(sources);
	Rcpp::DataFrame dfp = Rcpp::DataFrame(samples);

	// Check input data
	if(dfs.size() != 2 + (dfp.size() - 1) * 2)
	{
		stop("Wrong number of columns in source data frame");
	}

	// Check input data
	CharacterVector col_names = dfs.attr("names");
	CharacterVector col_names2 = dfp.attr("names");
	for( int i = 0 ; i < dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i];
		std::string name2 = "" + col_names2[i];
		if(name.compare(name2) != 0)
		{
			stop("Data frames must share column names");
		}
	}

	// Check input data
	for( int i = 1 ; i < dfp.size() ; i++ )
	{
		std::string name = "D" + col_names[i];
		std::string name2 = "" + col_names[i + dfp.size() - 1];
		if(name.compare(name2) != 0)
		{
			stop("Column " + name2 + " must be " + name + " in source data frame");
		}
	}

	// Check input data
	for( int i = dfp.size() ; i <= dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i + dfp.size() - 1];
		if (!(name.compare("n") == 0 || name.compare("N") == 0 || name.compare("num") == 0 || name.compare("Num") == 0))
		{
			stop("Last column in source data frame must be the number of samples");
		}
	}

	// Check input data
	std::string id = "" + col_names[0];
	if (!(id.compare("ID") == 0 || id.compare("Id") == 0 || id.compare("id") == 0))
	{
		stop("First data frame column must be ID");
	}

	int i, j, k, l, ii;
	int vars, nsource, points;
	double gof1, gof2, sum, avg, des, *max, *min, *try1, *try2;
	double **source, **source_d, **source_c, **point;
	int *source_n;

	// Number of variables
	vars = dfp.size() - 1;
	min = (double*) malloc( vars * sizeof(double) );
	max = (double*) malloc( vars * sizeof(double) );

	// Number of sources
	nsource = dfs.nrows();
	source   = (double**) malloc( nsource * sizeof(double*) );
	source_d = (double**) malloc( nsource * sizeof(double*) );
	source_c = (double**) malloc( nsource * sizeof(double*) );
	source_n = (int*) malloc( nsource * sizeof(int) );
	try1 = (double*) malloc( nsource * sizeof(double) );
	try2 = (double*) malloc( nsource * sizeof(double) );

	// Number of samples
	points = dfp.nrows();
	point = (double**) malloc( points * sizeof(double*) );

	// Output data frame
	//Rcpp::List myList(nsource * 2 + 2);
	Rcpp::List myList(nsource + 2);
	Rcpp::CharacterVector namevec;
	std::string name = "" + col_names[0];
	Rcpp::CharacterVector b(points*iter);// = dfp[name];
	myList[0] = b;
	namevec.push_back(name);
	Rcpp::CharacterVector c(points*iter);// = dfp[name];
	myList[1] = c;
	namevec.push_back("GOF");
	for(i = 0 ; i < nsource ; i++ )
	{
		std::string name = "" + col_names[1];
		Rcpp::CharacterVector b(points*iter);// = dfp[name];
		myList[i+2] = b;
		std::string name2 = "" + col_names[0];
		Rcpp::CharacterVector a = dfs[name2];
		namevec.push_back("w "+a[i]);

//		std::string name3 = "" + col_names[1];
//		Rcpp::CharacterVector c = dfp[name3];
//		myList[i+2+nsource] = c;
//		std::string name4 = "" + col_names[0];
//		Rcpp::CharacterVector d = dfs[name4];
//		namevec.push_back("Dw "+d[i]);
	}
	myList.attr("names") = namevec;

	// Read sources variables
	for(i = 0 ; i < nsource ; i++ )
	{
		source[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfs[name];
				source[i][j] = a[i];
		}
		source_d[i] = (double*) malloc( vars * sizeof(double) );
		source_c[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_d[i][j] = a[i];
		}
		for(j = vars ; j <= vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_n[i] = a[i];
		}
	}

	// Read samples variables
	for(i = 0 ; i < points ; i++ )
	{
		point[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfp[name];
				point[i][j] = a[i];
		}
	}

	// Struct to store trials
	type_try tried, best;
	tried.w = (double*) malloc( nsource * sizeof(double) );
	best.w = (double*) malloc( nsource * sizeof(double) );
	best.gof1 = 0.0;
	best.gof2 = 0.0;
	
	// central solution
	double *cw = (double*) malloc( nsource * sizeof(double) );
	// corrected mixture
	double *cm = (double*) malloc( vars * sizeof(double) );

	// Compute maximum values to normalize
	for(i = 0 ; i < vars ; i++ )
	{
		min[i] = max[i] = point[0][i];
		for(j = 1 ; j < points ; j++ )
		{
			if( point[j][i] < min[i] )
			{
				min[i] = point[j][i];
			}
			if( point[j][i] > max[i] )
			{
				max[i] = point[j][i];
			}
		}
		for(j = 0 ; j < nsource ; j++ )
		{
			if( source[j][i] < min[i] )
			{
				min[i] = source[j][i];
			}
			if( source[j][i] > max[i] )
			{
				max[i] = source[j][i];
			}
		}
	}

	// Init random seed
	rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	// test solutions iterations
	j = trials * 2;
	for( i = 2 ; i < nsource ; i++ )
	{
		j = j * trials;
	}
	trials = j;

	Progress p(points*iter, true);

	// Analize each point
	for( i = 0 ; i < points ; i++ )
	{
		// repeat iter times to get probability distribution
		for( ii = -1 ; ii < iter ; ii++ )
		{
			// screen
//			printf("\r                                            \r");
//			printf("> Analyzing %d of %d samples... %.1f%c",
//					i + 1, points, 100.0*(float)(ii+1)/(float)iter, '%');
//			fflush(stdout);

			if (Progress::check_abort() )
			return -1.0;

			p.increment();
		
			if(ii == -1 || ii == 1 || ii == 2)
			{
				// keep original mixture
				for( k = 0 ; k < vars ; k++ )
				{
					cm[k] = point[i][k];
				}
			}
			else
			{
				// compute corrected sources
				for( k = 0 ; k < vars ; k++ )
				{
					for( l = 0 ; l < nsource ; l++ )
					{
						source_c[l][k] = correct(source[l][k], source_d[l][k], source_n[l]);
					}
				}
				// compute corrected mixture
				for( k = 0 ; k < vars ; k++ )
				{
					sum = 0.0;
					for( l = 0 ; l < nsource ; l++ )
					{
						sum = sum + source_c[l][k] * cw[l];
					}
					cm[k] = sum;
				}
			}

			// test solutions
			for( j = 0 ; j < trials ; j++ )
			{
				// build random numbers wich sum 1
				// http://stats.stackexchange.com/questions/14059/generate-uniformly-distributed-weights-that-sum-to-unity
				for( k = 0 ; k < nsource - 1 ; k++ )
				{
					try2[k] = gsl_rng_uniform(rng);
				}
				qsort (try2, nsource - 1, sizeof(double), compare2);
				try1[0] = try2[0];
				for( k = 1 ; k < nsource - 1 ; k++ )
				{
					try1[k] = try2[k] - try2[k-1];
				}
				try1[nsource - 1] = 1.0 - try2[nsource - 2];

				// explore unphysical solutions to then discard
				for( k = 0 ; k < nsource ; k++ )
				{
					//try1[k] = try1[k] * 2.0 - 1.0 / (double)nsource;
					//try1[k] = try1[k] * 4.0 - 3.0 / (double)nsource; //increasing explored area
					//try1[k] = try1[k] * 6.0 - 5.0 / (double)nsource; //increasing explored area
					try1[k] = try1[k] * 10.0 - 9.0 / (double)nsource; //increasing explored area
				}

				// measure error
				gof1 = 0.0;
				gof2 = 0.0;
				for( k = 0 ; k < vars ; k++ )
				{
					sum = 0.0;
					for( l = 0 ; l < nsource ; l++ )
					{
						// Non corrected mean
						sum = sum + source[l][k] * try1[l];

						// Corrected mean
						//sum = sum + source_c[l][k] * try1[l];
					}
					// GOF1
					gof1 = gof1 + fabs( cm[k] - sum ) / ( max[k] - min[k] );

					// GOF2
					gof2 = gof2 + ( cm[k] - sum ) * ( cm[k] - sum ) / ( max[k] - min[k] ) / ( max[k] - min[k] );
				}
				gof1 = 1.0 - gof1 / (double) vars;
				gof2 = 1.0 - gof2 / (double) vars;

				// save best result
				if(gof2 > best.gof2 || j==0)
				{
					best.gof1 = gof1;
					best.gof2 = gof2;
					for( k = 0 ; k < nsource ; k++ )
					{
						best.w[k] = try1[k];
					}
				}

			}

			if(ii == -1)
			{
				// store central solution
				for( j = 0 ; j < nsource ; j++ )
				{
					cw[j] = best.w[j];
				}
			}
			else
			{
				// print average and standar deviation values in the top
				as<CharacterVector>(myList[0])[i*iter+ii] = i + 1;
				as<CharacterVector>(myList[1])[i*iter+ii] = best.gof1;
				for( j = 0 ; j < nsource ; j++ )
				{
					as<CharacterVector>(myList[2+j])[i*iter+ii] = best.w[j];
				}
			}
		}
	}
//	printf("\r                                            \r");

	Rcpp::DataFrame dfout(myList);
	return dfout;
}


//' least_squares_c
//'
//' @param sources Data frame containing sediment source samples
//' @param samples Data frame containing mixture samples
//' @param iter Iterations in the source variability analysis
//' @param seed Seed for the random number generator
//' @return Data frame containing the relative contribution solved by the least squares method
// [[Rcpp::export]]
Rcpp::DataFrame least_squares_c(SEXP sources, SEXP samples, int iter=100, int seed=69512)
{
	// Construct the data.frame object
	Rcpp::DataFrame dfs = Rcpp::DataFrame(sources);
	Rcpp::DataFrame dfp = Rcpp::DataFrame(samples);

	// Check input data
	if(dfs.size() != 2 + (dfp.size() - 1) * 2)
	{
		stop("Wrong number of columns in source data frame");
	}

	// Check input data
	CharacterVector col_names = dfs.attr("names");
	CharacterVector col_names2 = dfp.attr("names");
	for( int i = 0 ; i < dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i];
		std::string name2 = "" + col_names2[i];
		if(name.compare(name2) != 0)
		{
			stop("Data frames must share column names");
		}
	}

	// Check input data
	for( int i = 1 ; i < dfp.size() ; i++ )
	{
		std::string name = "D" + col_names[i];
		std::string name2 = "" + col_names[i + dfp.size() - 1];
		if(name.compare(name2) != 0)
		{
			stop("Column " + name2 + " must be " + name + " in source data frame");
		}
	}

	// Check input data
	for( int i = dfp.size() ; i <= dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i + dfp.size() - 1];
		if (!(name.compare("n") == 0 || name.compare("N") == 0 || name.compare("num") == 0 || name.compare("Num") == 0))
		{
			stop("Last column in source data frame must be the number of samples");
		}
	}

	// Check input data
	std::string id = "" + col_names[0];
	if (!(id.compare("ID") == 0 || id.compare("Id") == 0 || id.compare("id") == 0))
	{
		stop("First data frame column must be ID");
	}

	int i, j, k, l, ii;
	int vars, nsource, points;
	double gof1, gof2, sum, avg, des, *max, *min, *try1, *try2;
	double **source, **source_d, **source_c, **point;
	int *source_n;

	double min0, max0;
	double m1, m2, m3, m4;
	double s11, s21, s31, s41;
	double s12, s22, s32, s42;
	double s13, s23, s33, s43;
	double s14, s24, s34, s44;
	double det, w1, w2, w3, w4;
	double m11, m12, m13;
	double m21, m22, m23;
	double m31, m32, m33;
	double n11, n12, n13;
	double n21, n22, n23;
	double n31, n32, n33;
	double v1, v2, v3, v4;
	double o1, o2, o3, o4;
	double err = 0;
	
	// Number of variables
	vars = dfp.size() - 1;
	min = (double*) malloc( vars * sizeof(double) );
	max = (double*) malloc( vars * sizeof(double) );

	// Number of sources
	nsource = dfs.nrows();
	source   = (double**) malloc( nsource * sizeof(double*) );
	source_d = (double**) malloc( nsource * sizeof(double*) );
	source_c = (double**) malloc( nsource * sizeof(double*) );
	source_n = (int*) malloc( nsource * sizeof(int) );
	try1 = (double*) malloc( nsource * sizeof(double) );
	try2 = (double*) malloc( nsource * sizeof(double) );

	// Number of samples
	points = dfp.nrows();
	point = (double**) malloc( points * sizeof(double*) );

	// Check nsource and vars
	if (nsource == 3 && vars == 3)
	{
	}
	else if (nsource == 4 && vars == 4)
	{
	}
	else
	{
		stop("Least squares requires 3-3 or 4-4 tracers-sources");
	}

	// Output data frame
	//Rcpp::List myList(nsource * 2 + 2);
	Rcpp::List myList(nsource + 2);
	Rcpp::CharacterVector namevec;
	std::string name = "" + col_names[0];
	Rcpp::CharacterVector b(points*iter);// = dfp[name];
	myList[0] = b;
	namevec.push_back(name);
	Rcpp::CharacterVector c(points*iter);// = dfp[name];
	myList[1] = c;
	namevec.push_back("GOF");
	for(i = 0 ; i < nsource ; i++ )
	{
		std::string name = "" + col_names[1];
		Rcpp::CharacterVector b(points*iter);// = dfp[name];
		myList[i+2] = b;
		std::string name2 = "" + col_names[0];
		Rcpp::CharacterVector a = dfs[name2];
		namevec.push_back("w "+a[i]);

//		std::string name3 = "" + col_names[1];
//		Rcpp::CharacterVector c = dfp[name3];
//		myList[i+2+nsource] = c;
//		std::string name4 = "" + col_names[0];
//		Rcpp::CharacterVector d = dfs[name4];
//		namevec.push_back("Dw "+d[i]);
	}
	myList.attr("names") = namevec;

	// Read sources variables
	for(i = 0 ; i < nsource ; i++ )
	{
		source[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfs[name];
				source[i][j] = a[i];
		}
		source_d[i] = (double*) malloc( vars * sizeof(double) );
		source_c[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_d[i][j] = a[i];
		}
		for(j = vars ; j <= vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_n[i] = a[i];
		}
	}

	// Read samples variables
	for(i = 0 ; i < points ; i++ )
	{
		point[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfp[name];
				point[i][j] = a[i];
		}
	}

	// Struct to store trials
	type_try tried, best;
	tried.w = (double*) malloc( nsource * sizeof(double) );
	best.w = (double*) malloc( nsource * sizeof(double) );
	best.gof1 = 0.0;
	best.gof2 = 0.0;
	// Compute maximum values to normalize
	for(i = 0 ; i < vars ; i++ )
	{
		min[i] = 0.0;
		max[i] = 1.0;
	}

	// Init random seed
	rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	// test solutions iterations
//	j = trials;
//	for( i = 2 ; i < nsource ; i++ )
//	{
//		j = j * trials;
//	}
//	trials = j;

	//Progress p(points*iter, true);

	// Analize each point
	for( i = 0 ; i < points ; i++ )
	{
		// repeat iter times to get probability distribution
		for( ii = 0 ; ii < iter ; ii++ )
		{
			// screen
//			printf("\r                                            \r");
//			printf("> Analyzing %d of %d samples... %.1f%c",
//					i + 1, points, 100.0*(float)(ii+1)/(float)iter, '%');
//			fflush(stdout);

			//if (Progress::check_abort() )
			//return -1.0;

			//p.increment();
		
			// compute corrected sources
			for( k = 0 ; k < vars ; k++ )
			{
				for( l = 0 ; l < nsource ; l++ )
				{
					source_c[l][k] = correct(source[l][k], source_d[l][k], source_n[l]);
				}
			}

			// test solutions
			if(nsource == 3)
			{
				// tracer 1
				s11 = source_c[0][0];
				s12 = source_c[1][0];
				s13 = source_c[2][0];
				m1  = point[i][0];
				
				// tracer 2
				s21 = source_c[0][1];
				s22 = source_c[1][1];
				s23 = source_c[2][1];
				m2  = point[i][1];
				
				// tracer 3
				s31 = source_c[0][2];
				s32 = source_c[1][2];
				s33 = source_c[2][2];
				m3  = point[i][2];

				//printf("t1: %.4f %.4f %.4f\n", s11, s12, s13);
				//printf("t2: %.4f %.4f %.4f\n", s21, s22, s23);
				//printf("t3: %.4f %.4f %.4f\n", s31, s32, s33);
				//printf("\n");
				//printf("m:  %.4f %.4f %.4f\n", m1, m2, m3);
				//printf("\n");

				m11 = s33*s33-2.0*s31*s33+s31*s31+s23*s23-2.0*s21*s23+s21*s21+s13*s13-2.0*s11*s13+s11*s11;
				m12 = s33*s33-s32*s33-s31*s33+s31*s32+s23*s23-s22*s23-s21*s23+s21*s22+s13*s13-s12*s13-s11*s13+s11*s12;
				m21 = s33*s33-s32*s33-s31*s33+s31*s32+s23*s23-s22*s23-s21*s23+s21*s22+s13*s13-s12*s13-s11*s13+s11*s12;
				m22 = s33*s33-2.0*s32*s33+s32*s32+s23*s23-2.0*s22*s23+s22*s22+s13*s13-2.0*s12*s13+s12*s12;

				v1 = -(s33*m3-s31*m3+s23*m2-s21*m2+s13*m1-s11*m1-s33*s33+s31*s33-s23*s23+s21*s23-s13*s13+s11*s13);
				v2 = -(s33*m3-s32*m3+s23*m2-s22*m2+s13*m1-s12*m1-s33*s33+s32*s33-s23*s23+s22*s23-s13*s13+s12*s13);

				det = m11*m22-m12*m21;

				if(det != 0.0)
				{
					n11 = ( m22) / det;
					n12 = (-m12) / det;
					n21 = (-m21) / det;
					n22 = ( m11) / det;

					best.w[0] = n11 * v1 + n12 * v2;
					best.w[1] = n21 * v1 + n22 * v2;
					best.w[2] = 1.0 - best.w[0] - best.w[1];

					//printf("w:  %.4f %.4f %.4f\n", best.w[0], best.w[1], best.w[2]);

					o1 = s11 * best.w[0] + s12 * best.w[1] + s13 * best.w[2];
					o2 = s21 * best.w[0] + s22 * best.w[1] + s23 * best.w[2];
					o3 = s31 * best.w[0] + s32 * best.w[1] + s33 * best.w[2];

					err = (o1-m1)*(o1-m1)+(o2-m2)*(o2-m2)+(o3-m3)*(o3-m3);
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
					
					err = 1.E15;
				}
			}
			// test solutions
			if(nsource == 4)
			{
				// tracer 1
				s11 = source_c[0][0];
				s12 = source_c[1][0];
				s13 = source_c[2][0];
				s14 = source_c[3][0];
				m1  = point[i][0];
				
				// tracer 2
				s21 = source_c[0][1];
				s22 = source_c[1][1];
				s23 = source_c[2][1];
				s24 = source_c[3][1];
				m2  = point[i][1];
				
				// tracer 3
				s31 = source_c[0][2];
				s32 = source_c[1][2];
				s33 = source_c[2][2];
				s34 = source_c[3][2];
				m3  = point[i][2];
				
				// tracer 4
				s41 = source_c[0][3];
				s42 = source_c[1][3];
				s43 = source_c[2][3];
				s44 = source_c[3][3];
				m4  = point[i][3];

				m11 = s44*s44-2.0*s41*s44+s41*s41+s34*s34-2.0*s31*s34+s31*s31+s24*s24-2.0*s21*s24+s21*s21+s14*s14-2.0*s11*s14+s11*s11;
				m12 = s44*s44-s42*s44-s41*s44+s41*s42+s34*s34-s32*s34-s31*s34+s31*s32+s24*s24-s22*s24-s21*s24+s21*s22+s14*s14-s12*s14-s11*s14+s11*s12;
				m13 = s44*s44-s43*s44-s41*s44+s41*s43+s34*s34-s33*s34-s31*s34+s31*s33+s24*s24-s23*s24-s21*s24+s21*s23+s14*s14-s13*s14-s11*s14+s11*s13;
				m21 = s44*s44-s42*s44-s41*s44+s41*s42+s34*s34-s32*s34-s31*s34+s31*s32+s24*s24-s22*s24-s21*s24+s21*s22+s14*s14-s12*s14-s11*s14+s11*s12;
				m22 = s44*s44-2.0*s42*s44+s42*s42+s34*s34-2.0*s32*s34+s32*s32+s24*s24-2.0*s22*s24+s22*s22+s14*s14-2.0*s12*s14+s12*s12;
				m23 = s44*s44-s43*s44-s42*s44+s42*s43+s34*s34-s33*s34-s32*s34+s32*s33+s24*s24-s23*s24-s22*s24+s22*s23+s14*s14-s13*s14-s12*s14+s12*s13;
				m31 = s44*s44-s43*s44-s41*s44+s41*s43+s34*s34-s33*s34-s31*s34+s31*s33+s24*s24-s23*s24-s21*s24+s21*s23+s14*s14-s13*s14-s11*s14+s11*s13;
				m32 = s44*s44-s43*s44-s42*s44+s42*s43+s34*s34-s33*s34-s32*s34+s32*s33+s24*s24-s23*s24-s22*s24+s22*s23+s14*s14-s13*s14-s12*s14+s12*s13;
				m33 = s44*s44-2.0*s43*s44+s43*s43+s34*s34-2.0*s33*s34+s33*s33+s24*s24-2.0*s23*s24+s23*s23+s14*s14-2.0*s13*s14+s13*s13;

				v1 = -(s44*m4-s41*m4+s34*m3-s31*m3+s24*m2-s21*m2+s14*m1-s11*m1-s44*s44+s41*s44-s34*s34+s31*s34-s24*s24+s21*s24-s14*s14+s11*s14);
				v2 = -(s44*m4-s42*m4+s34*m3-s32*m3+s24*m2-s22*m2+s14*m1-s12*m1-s44*s44+s42*s44-s34*s34+s32*s34-s24*s24+s22*s24-s14*s14+s12*s14);
				v3 = -(s44*m4-s43*m4+s34*m3-s33*m3+s24*m2-s23*m2+s14*m1-s13*m1-s44*s44+s43*s44-s34*s34+s33*s34-s24*s24+s23*s24-s14*s14+s13*s14);

				det = (m11*m22-m12*m21)*m33+(m13*m21-m11*m23)*m32+(m12*m23-m13*m22)*m31;

				if(det != 0.0)
				{
					n11 = (m22*m33-m23*m32)/det;
					n12 = (m13*m32-m12*m33)/det;
					n13 = (m12*m23-m13*m22)/det;
					n21 = (m23*m31-m21*m33)/det;
					n22 = (m11*m33-m13*m31)/det;
					n23 = (m13*m21-m11*m23)/det;
					n31 = (m21*m32-m22*m31)/det;
					n32 = (m12*m31-m11*m32)/det;
					n33 = (m11*m22-m12*m21)/det;

					best.w[0] = n11 * v1 + n12 * v2 + n13 * v3;
					best.w[1] = n21 * v1 + n22 * v2 + n23 * v3;
					best.w[2] = n31 * v1 + n32 * v2 + n33 * v3;
					best.w[3] = 1.0 - best.w[0] - best.w[1] - best.w[2];

					o1 = s11 * best.w[0] + s12 * best.w[1] + s13 * best.w[2] + s14 * best.w[3];
					o2 = s21 * best.w[0] + s22 * best.w[1] + s23 * best.w[2] + s24 * best.w[3];
					o3 = s31 * best.w[0] + s32 * best.w[1] + s33 * best.w[2] + s34 * best.w[3];
					o4 = s41 * best.w[0] + s42 * best.w[1] + s43 * best.w[2] + s44 * best.w[3];

					err = (o1-m1)*(o1-m1)+(o2-m2)*(o2-m2)+(o3-m3)*(o3-m3)+(o4-m4)*(o4-m4);
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
					best.w[3] = 0.0;
					
					err = 1.E15;
				}
			}

			// print average and standar deviation values in the top
			as<CharacterVector>(myList[0])[i*iter+ii] = i + 1;
			as<CharacterVector>(myList[1])[i*iter+ii] = err;
			for( j = 0 ; j < nsource ; j++ )
			{
				as<CharacterVector>(myList[2+j])[i*iter+ii] = best.w[j];
			}
		}
	}
//	printf("\r                                            \r");

	Rcpp::DataFrame dfout(myList);
	return dfout;
}

//' triangles_random_c
//'
//' @param sources Data frame containing sediment source samples
//' @param samples Data frame containing mixture samples
//' @param tracer Tracer in which implement the function
//' @param iter Iterations in the source variability analysis
//' @param seed Seed for the random number generator
//' @return List of data frames containing all the possible prediction for each tracer
// [[Rcpp::export]]
Rcpp::DataFrame triangles_random_c(SEXP sources, SEXP samples, int tracer=0, int iter=100, int seed=69512)
{
	// Construct the data.frame object
	Rcpp::DataFrame dfs = Rcpp::DataFrame(sources);
	Rcpp::DataFrame dfp = Rcpp::DataFrame(samples);

	// Check input data
	if(dfs.size() != 2 + (dfp.size() - 1) * 2)
	{
		stop("Wrong number of columns in source data frame");
	}

	// Check input data
	CharacterVector col_names = dfs.attr("names");
	CharacterVector col_names2 = dfp.attr("names");
	for( int i = 0 ; i < dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i];
		std::string name2 = "" + col_names2[i];
		if(name.compare(name2) != 0)
		{
			stop("Data frames must share column names");
		}
	}

	// Check input data
	for( int i = 1 ; i < dfp.size() ; i++ )
	{
		std::string name = "D" + col_names[i];
		std::string name2 = "" + col_names[i + dfp.size() - 1];
		if(name.compare(name2) != 0)
		{
			stop("Column " + name2 + " must be " + name + " in source data frame");
		}
	}

	// Check input data
	for( int i = dfp.size() ; i <= dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i + dfp.size() - 1];
		if (!(name.compare("n") == 0 || name.compare("N") == 0 || name.compare("num") == 0 || name.compare("Num") == 0))
		{
			stop("Last column in source data frame must be the number of samples");
		}
	}

	// Check input data
	std::string id = "" + col_names[0];
	if (!(id.compare("ID") == 0 || id.compare("Id") == 0 || id.compare("id") == 0))
	{
		stop("First data frame column must be ID");
	}

	int i, j, k, l, ii;
	int vars, nsource, points;
	double gof1, gof2, sum, avg, des, *max, *min, *try1, *try2;
	double **source, **source_d, **source_c, **point;
	int *source_n;

	double min0, max0;
	double m1, m2, m3;
	double s11, s21, s31, s41;
	double s12, s22, s32, s42;
	double s13, s23, s33, s43;
	double det, w1, w2, w3, w4;

	// Number of variables
	vars = dfp.size() - 1;
	min = (double*) malloc( vars * sizeof(double) );
	max = (double*) malloc( vars * sizeof(double) );

	// Number of sources
	nsource = dfs.nrows();
	source   = (double**) malloc( nsource * sizeof(double*) );
	source_d = (double**) malloc( nsource * sizeof(double*) );
	source_c = (double**) malloc( nsource * sizeof(double*) );
	source_n = (int*) malloc( nsource * sizeof(int) );
	try1 = (double*) malloc( nsource * sizeof(double) );
	try2 = (double*) malloc( nsource * sizeof(double) );

	// Number of samples
	points = dfp.nrows();
	point = (double**) malloc( points * sizeof(double*) );

	// Check tracer
	if (tracer < 0 || tracer >= vars)
	{
		stop("Requested tracer is out of range");
	}

	// Check nsource
	if (!(nsource == 3 || nsource == 4))
	{
		stop("Number of sources not implemented");
	}

	// Output data frame
	//Rcpp::List myList(nsource * 2 + 2);
	Rcpp::List myList(nsource + 2);
	Rcpp::CharacterVector namevec;
	std::string name = "" + col_names[0];
	Rcpp::CharacterVector b(points*iter);// = dfp[name];
	myList[0] = b;
	namevec.push_back(name);
	Rcpp::CharacterVector c(points*iter);// = dfp[name];
	myList[1] = c;
	namevec.push_back("GOF");
	for(i = 0 ; i < nsource ; i++ )
	{
		std::string name = "" + col_names[1];
		Rcpp::CharacterVector b(points*iter);// = dfp[name];
		myList[i+2] = b;
		std::string name2 = "" + col_names[0];
		Rcpp::CharacterVector a = dfs[name2];
		namevec.push_back("w "+a[i]);

//		std::string name3 = "" + col_names[1];
//		Rcpp::CharacterVector c = dfp[name3];
//		myList[i+2+nsource] = c;
//		std::string name4 = "" + col_names[0];
//		Rcpp::CharacterVector d = dfs[name4];
//		namevec.push_back("Dw "+d[i]);
	}
	myList.attr("names") = namevec;

	// Read sources variables
	for(i = 0 ; i < nsource ; i++ )
	{
		source[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfs[name];
				source[i][j] = a[i];
		}
		source_d[i] = (double*) malloc( vars * sizeof(double) );
		source_c[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_d[i][j] = a[i];
		}
		for(j = vars ; j <= vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_n[i] = a[i];
		}
	}

	// Read samples variables
	for(i = 0 ; i < points ; i++ )
	{
		point[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfp[name];
				point[i][j] = a[i];
		}
	}

	// Struct to store trials
	type_try tried, best;
	tried.w = (double*) malloc( nsource * sizeof(double) );
	best.w = (double*) malloc( nsource * sizeof(double) );
	best.gof1 = 0.0;
	best.gof2 = 0.0;
	// Compute maximum values to normalize
	for(i = 0 ; i < vars ; i++ )
	{
		min[i] = max[i] = point[0][i];
		for(j = 1 ; j < points ; j++ )
		{
			if( point[j][i] < min[i] )
			{
				min[i] = point[j][i];
			}
			if( point[j][i] > max[i] )
			{
				max[i] = point[j][i];
			}
		}
		for(j = 0 ; j < nsource ; j++ )
		{
			if( source[j][i] < min[i] )
			{
				min[i] = source[j][i];
			}
			if( source[j][i] > max[i] )
			{
				max[i] = source[j][i];
			}
		}
	}

	// Init random seed
	rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	// test solutions iterations
//	j = trials;
//	for( i = 2 ; i < nsource ; i++ )
//	{
//		j = j * trials;
//	}
//	trials = j;

	Progress p(points*iter, true);

	// Analize each point
	for( i = 0 ; i < points ; i++ )
	{
		// repeat iter times to get probability distribution
		for( ii = 0 ; ii < iter ; ii++ )
		{
			// screen
//			printf("\r                                            \r");
//			printf("> Analyzing %d of %d samples... %.1f%c",
//					i + 1, points, 100.0*(float)(ii+1)/(float)iter, '%');
//			fflush(stdout);

			if (Progress::check_abort() )
			return -1.0;

			p.increment();
		
			// compute corrected sources
			for( k = 0 ; k < vars ; k++ )
			{
				for( l = 0 ; l < nsource ; l++ )
				{
					source_c[l][k] = correct(source[l][k], source_d[l][k], source_n[l]);
				}
			}

			// test solutions
			if(nsource == 3)
			{
				// requested tracer
				s11 = source_c[0][tracer];
				s21 = source_c[1][tracer];
				s31 = source_c[2][tracer];
				m1  = point[i][tracer];

				// select random tracer
				k = tracer;
				while(k == tracer)
				{
					k = (int)(gsl_rng_uniform(rng) * (double)vars);
				}
				s12 = source_c[0][k];
				s22 = source_c[1][k];
				s32 = source_c[2][k];
				m2  = point[i][k];

				// build random tracer
//				max0 = 0.0 ; min0 = 1.0;
//				s12 = gsl_rng_uniform(rng); if(s12 > max0){ max0 = s12; } if(s12 < min0){ min0 = s12; }
//				s22 = gsl_rng_uniform(rng); if(s22 > max0){ max0 = s22; } if(s22 < min0){ min0 = s22; }
//				s32 = gsl_rng_uniform(rng); if(s32 > max0){ max0 = s32; } if(s32 < min0){ min0 = s32; }
//				m2 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				det = -s22*s31+s12*s31+s21*s32-s11*s32+s11*s22-s12*s21;

				if(det != 0.0)
				{
					w1 = s21*s32-s22*s31 + m1*(s22-s32) + m2*(s31-s21);
					w2 = s12*s31-s11*s32 + m1*(s32-s12) + m2*(s11-s31);
					w3 = s11*s22-s12*s21 + m1*(s12-s22) + m2*(s21-s11);

					best.w[0] = w1 / det;
					best.w[1] = w2 / det;
					best.w[2] = w3 / det;
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
				}
			}
			// test solutions
			if(nsource == 4)
			{
				// requested tracer
				s11 = source_c[0][tracer];
				s21 = source_c[1][tracer];
				s31 = source_c[2][tracer];
				s41 = source_c[3][tracer];
				m1  = point[i][tracer];

				// select random tracer
				k = tracer;
				while(k == tracer)
				{
					k = (int)(gsl_rng_uniform(rng) * (double)vars);
				}
				s12 = source_c[0][k];
				s22 = source_c[1][k];
				s32 = source_c[2][k];
				s42 = source_c[3][k];
				m2  = point[i][k];

				// build random tracer
//				max0 = 0.0 ; min0 = 1.0;
//				s12 = gsl_rng_uniform(rng); if(s12 > max0){ max0 = s12; } if(s12 < min0){ min0 = s12; }
//				s22 = gsl_rng_uniform(rng); if(s22 > max0){ max0 = s22; } if(s22 < min0){ min0 = s22; }
//				s32 = gsl_rng_uniform(rng); if(s32 > max0){ max0 = s32; } if(s32 < min0){ min0 = s32; }
//				s42 = gsl_rng_uniform(rng); if(s42 > max0){ max0 = s42; } if(s42 < min0){ min0 = s42; }
//				m2 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				// select random tracer
				l = tracer;
				while(l == tracer || l == k)
				{
					l = (int)(gsl_rng_uniform(rng) * (double)vars);
				}
				s13 = source_c[0][l];
				s23 = source_c[1][l];
				s33 = source_c[2][l];
				s43 = source_c[3][l];
				m3  = point[i][l];

				// build random tracer
//				max0 = 0.0 ; min0 = 1.0;
//				s13 = gsl_rng_uniform(rng); if(s13 > max0){ max0 = s13; } if(s13 < min0){ min0 = s13; }
//				s23 = gsl_rng_uniform(rng); if(s23 > max0){ max0 = s23; } if(s23 < min0){ min0 = s23; }
//				s33 = gsl_rng_uniform(rng); if(s33 > max0){ max0 = s33; } if(s33 < min0){ min0 = s33; }
//				s43 = gsl_rng_uniform(rng); if(s43 > max0){ max0 = s43; } if(s43 < min0){ min0 = s43; }
//				m3 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				det = s21*(s32*s43-s33*s42)-s11*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+s11*(s22*s43-s23*s42)+s31*(s12*s43-s13*s42)-s21*(s12*s43-s13*s42)+(s22*s33-s23*s32)*s41-(s12*s33-s13*s32)*s41+(s12*s23-s13*s22)*s41-s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31;

				if(det != 0.0)
				{
					w1 = s21*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+(s22*s33-s23*s32)*s41 + m1*(-s32*s43+s22*s43+s33*s42-s23*s42-s22*s33+s23*s32) + m2*(s31*s43-s21*s43-s33*s41+s23*s41+s21*s33-s23*s31) + m3*(-s31*s42+s21*s42+s32*s41-s22*s41-s21*s32+s22*s31);

					w2 = -s11*(s32*s43-s33*s42)+s31*(s12*s43-s13*s42)-(s12*s33-s13*s32)*s41 + m1*(s32*s43-s12*s43-s33*s42+s13*s42+s12*s33-s13*s32) + m2*(-s31*s43+s11*s43+s33*s41-s13*s41-s11*s33+s13*s31) + m3*(s31*s42-s11*s42-s32*s41+s12*s41+s11*s32-s12*s31);

					w3 = s11*(s22*s43-s23*s42)-s21*(s12*s43-s13*s42)+(s12*s23-s13*s22)*s41 + m1*(-s22*s43+s12*s43+s23*s42-s13*s42-s12*s23+s13*s22) + m2*(s21*s43-s11*s43-s23*s41+s13*s41+s11*s23-s13*s21) + m3*(-s21*s42+s11*s42+s22*s41-s12*s41-s11*s22+s12*s21);

					w4 = -s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31 + m1*(s22*s33-s12*s33-s23*s32+s13*s32+s12*s23-s13*s22) + m2*(-s21*s33+s11*s33+s23*s31-s13*s31-s11*s23+s13*s21) + m3*(s21*s32-s11*s32-s22*s31+s12*s31+s11*s22-s12*s21);

					best.w[0] = w1 / det;
					best.w[1] = w2 / det;
					best.w[2] = w3 / det;
					best.w[3] = w4 / det;
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
					best.w[3] = 0.0;
				}
			}

			// print average and standar deviation values in the top
			as<CharacterVector>(myList[0])[i*iter+ii] = i + 1;
			as<CharacterVector>(myList[1])[i*iter+ii] = 1.0;
			for( j = 0 ; j < nsource ; j++ )
			{
				as<CharacterVector>(myList[2+j])[i*iter+ii] = best.w[j];
			}
		}
	}
//	printf("\r                                            \r");

	Rcpp::DataFrame dfout(myList);
	return dfout;
}

//' triangles_virtual_c
//'
//' @param sources Data frame containing sediment source samples
//' @param samples Data frame containing mixture samples
//' @param tracer Tracer in which implement the function
//' @param iter Iterations in the source variability analysis
//' @param seed Seed for the random number generator
//' @return List of data frames containing all the possible prediction for each tracer inside the dataset
// [[Rcpp::export]]
Rcpp::DataFrame triangles_virtual_c(SEXP sources, SEXP samples, int tracer=0, int iter=100, int seed=69512)
{
	// Construct the data.frame object
	Rcpp::DataFrame dfs = Rcpp::DataFrame(sources);
	Rcpp::DataFrame dfp = Rcpp::DataFrame(samples);

	// Check input data
	if(dfs.size() != 2 + (dfp.size() - 1) * 2)
	{
		stop("Wrong number of columns in source data frame");
	}

	// Check input data
	CharacterVector col_names = dfs.attr("names");
	CharacterVector col_names2 = dfp.attr("names");
	for( int i = 0 ; i < dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i];
		std::string name2 = "" + col_names2[i];
		if(name.compare(name2) != 0)
		{
			stop("Data frames must share column names");
		}
	}

	// Check input data
	for( int i = 1 ; i < dfp.size() ; i++ )
	{
		std::string name = "D" + col_names[i];
		std::string name2 = "" + col_names[i + dfp.size() - 1];
		if(name.compare(name2) != 0)
		{
			stop("Column " + name2 + " must be " + name + " in source data frame");
		}
	}

	// Check input data
	for( int i = dfp.size() ; i <= dfp.size() ; i++ )
	{
		std::string name = "" + col_names[i + dfp.size() - 1];
		if (!(name.compare("n") == 0 || name.compare("N") == 0 || name.compare("num") == 0 || name.compare("Num") == 0))
		{
			stop("Last column in source data frame must be the number of samples");
		}
	}

	// Check input data
	std::string id = "" + col_names[0];
	if (!(id.compare("ID") == 0 || id.compare("Id") == 0 || id.compare("id") == 0))
	{
		stop("First data frame column must be ID");
	}

	int i, j, k, l, ii;
	int vars, nsource, points;
	double gof1, gof2, sum, avg, des, *max, *min, *try1, *try2;
	double **source, **source_d, **source_c, **point;
	int *source_n;

	double min0, max0;
	double m1, m2, m3, m4;
	double s11, s21, s31, s41, s51;
	double s12, s22, s32, s42, s52;
	double s13, s23, s33, s43, s53;
	double s14, s24, s34, s44, s54;
	double det, w1, w2, w3, w4, w5;

	// Number of variables
	vars = dfp.size() - 1;
	min = (double*) malloc( vars * sizeof(double) );
	max = (double*) malloc( vars * sizeof(double) );

	// Number of sources
	nsource = dfs.nrows();
	source   = (double**) malloc( nsource * sizeof(double*) );
	source_d = (double**) malloc( nsource * sizeof(double*) );
	source_c = (double**) malloc( nsource * sizeof(double*) );
	source_n = (int*) malloc( nsource * sizeof(int) );
	try1 = (double*) malloc( nsource * sizeof(double) );
	try2 = (double*) malloc( nsource * sizeof(double) );

	// Number of samples
	points = dfp.nrows();
	point = (double**) malloc( points * sizeof(double*) );

	// Check tracer
	if (tracer < 0 || tracer >= vars)
	{
		stop("Requested tracer is out of range");
	}

	// Check nsource
	if (!(nsource == 2 || nsource == 3 || nsource == 4 || nsource == 5))
	{		
		stop("Number of sources not implemented");
	}

	// Output data frame
	//Rcpp::List myList(nsource * 2 + 2);
	Rcpp::List myList(nsource + 2);
	Rcpp::CharacterVector namevec;
	std::string name = "" + col_names[0];
	Rcpp::CharacterVector b(points*iter);// = dfp[name];
	myList[0] = b;
	namevec.push_back(name);
	Rcpp::CharacterVector c(points*iter);// = dfp[name];
	myList[1] = c;
	namevec.push_back("GOF");
	for(i = 0 ; i < nsource ; i++ )
	{
		std::string name = "" + col_names[1];
		Rcpp::CharacterVector b(points*iter);// = dfp[name];
		myList[i+2] = b;
		std::string name2 = "" + col_names[0];
		Rcpp::CharacterVector a = dfs[name2];
		namevec.push_back("w "+a[i]);

//		std::string name3 = "" + col_names[1];
//		Rcpp::CharacterVector c = dfp[name3];
//		myList[i+2+nsource] = c;
//		std::string name4 = "" + col_names[0];
//		Rcpp::CharacterVector d = dfs[name4];
//		namevec.push_back("Dw "+d[i]);
	}
	myList.attr("names") = namevec;

	// Read sources variables
	for(i = 0 ; i < nsource ; i++ )
	{
		source[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfs[name];
				source[i][j] = a[i];
		}
		source_d[i] = (double*) malloc( vars * sizeof(double) );
		source_c[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_d[i][j] = a[i];
		}
		for(j = vars ; j <= vars ; j++ )
		{
				std::string name = "" + col_names[j+1+vars];
				Rcpp::NumericVector a = dfs[name];
				source_n[i] = a[i];
		}
	}

	// Read samples variables
	for(i = 0 ; i < points ; i++ )
	{
		point[i] = (double*) malloc( vars * sizeof(double) );
		for(j = 0 ; j < vars ; j++ )
		{
				std::string name = "" + col_names[j+1];
				Rcpp::NumericVector a = dfp[name];
				point[i][j] = a[i];
		}
	}

	// Struct to store trials
	type_try tried, best;
	tried.w = (double*) malloc( nsource * sizeof(double) );
	best.w = (double*) malloc( nsource * sizeof(double) );
	best.gof1 = 0.0;
	best.gof2 = 0.0;
	// Compute maximum values to normalize
	for(i = 0 ; i < vars ; i++ )
	{
		min[i] = max[i] = point[0][i];
		for(j = 1 ; j < points ; j++ )
		{
			if( point[j][i] < min[i] )
			{
				min[i] = point[j][i];
			}
			if( point[j][i] > max[i] )
			{
				max[i] = point[j][i];
			}
		}
		for(j = 0 ; j < nsource ; j++ )
		{
			if( source[j][i] < min[i] )
			{
				min[i] = source[j][i];
			}
			if( source[j][i] > max[i] )
			{
				max[i] = source[j][i];
			}
		}
	}

	// Init random seed
	rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	// test solutions iterations
//	j = trials;
//	for( i = 2 ; i < nsource ; i++ )
//	{
//		j = j * trials;
//	}
//	trials = j;

	Progress p(points*iter, true);

	// Analize each point
	for( i = 0 ; i < points ; i++ )
	{
		// repeat iter times to get probability distribution
		for( ii = 0 ; ii < iter ; ii++ )
		{
			// screen
//			printf("\r                                            \r");
//			printf("> Analyzing %d of %d samples... %.1f%c",
//					i + 1, points, 100.0*(float)(ii+1)/(float)iter, '%');
//			fflush(stdout);

			if (Progress::check_abort() )
			return -1.0;

			p.increment();
		
			// compute corrected sources
			for( k = 0 ; k < vars ; k++ )
			{
				for( l = 0 ; l < nsource ; l++ )
				{
					source_c[l][k] = correct(source[l][k], source_d[l][k], source_n[l]);
				}
			}

			// test solutions	
			if(nsource == 2)	
			{	
				// requested tracer	
				s11 = source_c[0][tracer];	
				s21 = source_c[1][tracer];	
				m1  = point[i][tracer];	
				det = s11 - s21;	
				if(det != 0.0)	
				{	
					w1 = (m1 - s21) / det;	
					w2 = 1.0 - w1;	
					best.w[0] = w1;	
					best.w[1] = w2;	
				}	
				else	
				{	
					best.w[0] = 0.0;	
					best.w[1] = 0.0;	
				}	
			}
			// test solutions
			if(nsource == 3)
			{
				// requested tracer
				s11 = source_c[0][tracer];
				s21 = source_c[1][tracer];
				s31 = source_c[2][tracer];
				m1  = point[i][tracer];

				// select random tracer
//				k = tracer;
//				while(k == tracer)
//				{
//					k = (int)(gsl_rng_uniform(rng) * (double)vars);
//				}
//				s12 = source_c[0][k];
//				s22 = source_c[1][k];
//				s32 = source_c[2][k];
//				m2  = point[i][k];

				// build random tracer
				max0 = 0.0 ; min0 = 1.0;
				s12 = gsl_rng_uniform(rng); if(s12 > max0){ max0 = s12; } if(s12 < min0){ min0 = s12; }
				s22 = gsl_rng_uniform(rng); if(s22 > max0){ max0 = s22; } if(s22 < min0){ min0 = s22; }
				s32 = gsl_rng_uniform(rng); if(s32 > max0){ max0 = s32; } if(s32 < min0){ min0 = s32; }
				m2 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				det = -s22*s31+s12*s31+s21*s32-s11*s32+s11*s22-s12*s21;

				if(det != 0.0)
				{
					w1 = s21*s32-s22*s31 + m1*(s22-s32) + m2*(s31-s21);
					w2 = s12*s31-s11*s32 + m1*(s32-s12) + m2*(s11-s31);
					w3 = s11*s22-s12*s21 + m1*(s12-s22) + m2*(s21-s11);

					best.w[0] = w1 / det;
					best.w[1] = w2 / det;
					best.w[2] = w3 / det;
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
				}
			}
			// test solutions
			if(nsource == 4)
			{
				// requested tracer
				s11 = source_c[0][tracer];
				s21 = source_c[1][tracer];
				s31 = source_c[2][tracer];
				s41 = source_c[3][tracer];
				m1  = point[i][tracer];

				// select random tracer
//				k = tracer;
//				while(k == tracer)
//				{
//					k = (int)(gsl_rng_uniform(rng) * (double)vars);
//				}
//				s12 = source_c[0][k];
//				s22 = source_c[1][k];
//				s32 = source_c[2][k];
//				s42 = source_c[3][k];
//				m2  = point[i][k];

				// build random tracer
				max0 = 0.0 ; min0 = 1.0;
				s12 = gsl_rng_uniform(rng); if(s12 > max0){ max0 = s12; } if(s12 < min0){ min0 = s12; }
				s22 = gsl_rng_uniform(rng); if(s22 > max0){ max0 = s22; } if(s22 < min0){ min0 = s22; }
				s32 = gsl_rng_uniform(rng); if(s32 > max0){ max0 = s32; } if(s32 < min0){ min0 = s32; }
				s42 = gsl_rng_uniform(rng); if(s42 > max0){ max0 = s42; } if(s42 < min0){ min0 = s42; }
				m2 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				// select random tracer
//				l = tracer;
//				while(l == tracer || l == k)
//				{
//					k = (int)(gsl_rng_uniform(rng) * (double)vars);
//				}
//				s13 = source_c[0][l];
//				s23 = source_c[1][l];
//				s33 = source_c[2][l];
//				s43 = source_c[3][l];
//				m3  = point[i][l];

				// build random tracer
				max0 = 0.0 ; min0 = 1.0;
				s13 = gsl_rng_uniform(rng); if(s13 > max0){ max0 = s13; } if(s13 < min0){ min0 = s13; }
				s23 = gsl_rng_uniform(rng); if(s23 > max0){ max0 = s23; } if(s23 < min0){ min0 = s23; }
				s33 = gsl_rng_uniform(rng); if(s33 > max0){ max0 = s33; } if(s33 < min0){ min0 = s33; }
				s43 = gsl_rng_uniform(rng); if(s43 > max0){ max0 = s43; } if(s43 < min0){ min0 = s43; }
				m3 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				det = s21*(s32*s43-s33*s42)-s11*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+s11*(s22*s43-s23*s42)+s31*(s12*s43-s13*s42)-s21*(s12*s43-s13*s42)+(s22*s33-s23*s32)*s41-(s12*s33-s13*s32)*s41+(s12*s23-s13*s22)*s41-s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31;

				if(det != 0.0)
				{
					w1 = s21*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+(s22*s33-s23*s32)*s41 + m1*(-s32*s43+s22*s43+s33*s42-s23*s42-s22*s33+s23*s32) + m2*(s31*s43-s21*s43-s33*s41+s23*s41+s21*s33-s23*s31) + m3*(-s31*s42+s21*s42+s32*s41-s22*s41-s21*s32+s22*s31);

					w2 = -s11*(s32*s43-s33*s42)+s31*(s12*s43-s13*s42)-(s12*s33-s13*s32)*s41 + m1*(s32*s43-s12*s43-s33*s42+s13*s42+s12*s33-s13*s32) + m2*(-s31*s43+s11*s43+s33*s41-s13*s41-s11*s33+s13*s31) + m3*(s31*s42-s11*s42-s32*s41+s12*s41+s11*s32-s12*s31);

					w3 = s11*(s22*s43-s23*s42)-s21*(s12*s43-s13*s42)+(s12*s23-s13*s22)*s41 + m1*(-s22*s43+s12*s43+s23*s42-s13*s42-s12*s23+s13*s22) + m2*(s21*s43-s11*s43-s23*s41+s13*s41+s11*s23-s13*s21) + m3*(-s21*s42+s11*s42+s22*s41-s12*s41-s11*s22+s12*s21);

					w4 = -s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31 + m1*(s22*s33-s12*s33-s23*s32+s13*s32+s12*s23-s13*s22) + m2*(-s21*s33+s11*s33+s23*s31-s13*s31-s11*s23+s13*s21) + m3*(s21*s32-s11*s32-s22*s31+s12*s31+s11*s22-s12*s21);

					best.w[0] = w1 / det;
					best.w[1] = w2 / det;
					best.w[2] = w3 / det;
					best.w[3] = w4 / det;
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
					best.w[3] = 0.0;
				}
			}
			// test solutions
			if(nsource == 5)
			{
				// requested tracer
				s11 = source_c[0][tracer];
				s21 = source_c[1][tracer];
				s31 = source_c[2][tracer];
				s41 = source_c[3][tracer];
				s51 = source_c[4][tracer];
				m1  = point[i][tracer];

				// build random tracer
				max0 = 0.0 ; min0 = 1.0;
				s12 = gsl_rng_uniform(rng); if(s12 > max0){ max0 = s12; } if(s12 < min0){ min0 = s12; }
				s22 = gsl_rng_uniform(rng); if(s22 > max0){ max0 = s22; } if(s22 < min0){ min0 = s22; }
				s32 = gsl_rng_uniform(rng); if(s32 > max0){ max0 = s32; } if(s32 < min0){ min0 = s32; }
				s42 = gsl_rng_uniform(rng); if(s42 > max0){ max0 = s42; } if(s42 < min0){ min0 = s42; }
				s52 = gsl_rng_uniform(rng); if(s52 > max0){ max0 = s52; } if(s52 < min0){ min0 = s52; }
				m2 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				// build random tracer
				max0 = 0.0 ; min0 = 1.0;
				s13 = gsl_rng_uniform(rng); if(s13 > max0){ max0 = s13; } if(s13 < min0){ min0 = s13; }
				s23 = gsl_rng_uniform(rng); if(s23 > max0){ max0 = s23; } if(s23 < min0){ min0 = s23; }
				s33 = gsl_rng_uniform(rng); if(s33 > max0){ max0 = s33; } if(s33 < min0){ min0 = s33; }
				s43 = gsl_rng_uniform(rng); if(s43 > max0){ max0 = s43; } if(s43 < min0){ min0 = s43; }
				s53 = gsl_rng_uniform(rng); if(s53 > max0){ max0 = s53; } if(s53 < min0){ min0 = s53; }
				m3 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				// build random tracer
				max0 = 0.0 ; min0 = 1.0;
				s14 = gsl_rng_uniform(rng); if(s14 > max0){ max0 = s14; } if(s14 < min0){ min0 = s14; }
				s24 = gsl_rng_uniform(rng); if(s24 > max0){ max0 = s24; } if(s24 < min0){ min0 = s24; }
				s34 = gsl_rng_uniform(rng); if(s34 > max0){ max0 = s34; } if(s34 < min0){ min0 = s34; }
				s44 = gsl_rng_uniform(rng); if(s44 > max0){ max0 = s44; } if(s44 < min0){ min0 = s44; }
				s54 = gsl_rng_uniform(rng); if(s54 > max0){ max0 = s54; } if(s54 < min0){ min0 = s54; }
				m4 = min0 + gsl_rng_uniform(rng) * (max0 - min0);

				// maxima script
				// M:matrix([1,1,1,1,1],[s11,s21,s31,s41,s51],[s12,s22,s32,s42,s52],[s13,s23,s33,s43,s53],[s14,s24,s34,s44,s54]);
				// determinant(M);
				// ratsimp(invert(M)*determinant(M));

				det = s21*(s32*(s43*s54-s44*s53)-s42*(s33*s54-s34*s53)+(s33*s44-s34*s43)*s52)-s11*(s32*(s43*s54-s44*s53)-s42*(s33*s54-s34*s53)+(s33*s44-s34*s43)*s52)-s31*(s22*(s43*s54-s44*s53)-s42*(s23*s54-s24*s53)+(s23*s44-s24*s43)*s52)+s11*(s22*(s43*s54-s44*s53)-s42*(s23*s54-s24*s53)+(s23*s44-s24*s43)*s52)+s31*(s12*(s43*s54-s44*s53)-s42*(s13*s54-s14*s53)+(s13*s44-s14*s43)*s52)-s21*(s12*(s43*s54-s44*s53)-s42*(s13*s54-s14*s53)+(s13*s44-s14*s43)*s52)+s41*(s22*(s33*s54-s34*s53)-s32*(s23*s54-s24*s53)+(s23*s34-s24*s33)*s52)-s11*(s22*(s33*s54-s34*s53)-s32*(s23*s54-s24*s53)+(s23*s34-s24*s33)*s52)-s41*(s12*(s33*s54-s34*s53)-s32*(s13*s54-s14*s53)+(s13*s34-s14*s33)*s52)+s21*(s12*(s33*s54-s34*s53)-s32*(s13*s54-s14*s53)+(s13*s34-s14*s33)*s52)+s41*(s12*(s23*s54-s24*s53)-s22*(s13*s54-s14*s53)+(s13*s24-s14*s23)*s52)-s31*(s12*(s23*s54-s24*s53)-s22*(s13*s54-s14*s53)+(s13*s24-s14*s23)*s52)-(s22*(s33*s44-s34*s43)-s32*(s23*s44-s24*s43)+(s23*s34-s24*s33)*s42)*s51+(s12*(s33*s44-s34*s43)-s32*(s13*s44-s14*s43)+(s13*s34-s14*s33)*s42)*s51-(s12*(s23*s44-s24*s43)-s22*(s13*s44-s14*s43)+(s13*s24-s14*s23)*s42)*s51+(s12*(s23*s34-s24*s33)-s22*(s13*s34-s14*s33)+(s13*s24-s14*s23)*s32)*s51+s11*(s22*(s33*s44-s34*s43)-s32*(s23*s44-s24*s43)+(s23*s34-s24*s33)*s42)-s21*(s12*(s33*s44-s34*s43)-s32*(s13*s44-s14*s43)+(s13*s34-s14*s33)*s42)+s31*(s12*(s23*s44-s24*s43)-s22*(s13*s44-s14*s43)+(s13*s24-s14*s23)*s42)-(s12*(s23*s34-s24*s33)-s22*(s13*s34-s14*s33)+(s13*s24-s14*s23)*s32)*s41;
	
				if(det != 0.0)
				{
					w1 = m1*(((s22-s32)*s43+(s33-s23)*s42-s22*s33+s23*s32)*s54+((s32-s22)*s44+(s24-s34)*s42+s22*s34-s24*s32)*s53+((s23-s33)*s44+(s34-s24)*s43-s23*s34+s24*s33)*s52+(s22*s33-s23*s32)*s44+(s24*s32-s22*s34)*s43+(s23*s34-s24*s33)*s42)+m2*(((s31-s21)*s43+(s23-s33)*s41+s21*s33-s23*s31)*s54+((s21-s31)*s44+(s34-s24)*s41-s21*s34+s24*s31)*s53+((s33-s23)*s44+(s24-s34)*s43+s23*s34-s24*s33)*s51+(s23*s31-s21*s33)*s44+(s21*s34-s24*s31)*s43+(s24*s33-s23*s34)*s41)+m3*(((s21-s31)*s42+(s32-s22)*s41-s21*s32+s22*s31)*s54+((s31-s21)*s44+(s24-s34)*s41+s21*s34-s24*s31)*s52+((s22-s32)*s44+(s34-s24)*s42-s22*s34+s24*s32)*s51+(s21*s32-s22*s31)*s44+(s24*s31-s21*s34)*s42+(s22*s34-s24*s32)*s41)+((s21*s32-s22*s31)*s43+(s23*s31-s21*s33)*s42+(s22*s33-s23*s32)*s41)*s54+m4*(((s31-s21)*s42+(s22-s32)*s41+s21*s32-s22*s31)*s53+((s21-s31)*s43+(s33-s23)*s41-s21*s33+s23*s31)*s52+((s32-s22)*s43+(s23-s33)*s42+s22*s33-s23*s32)*s51+(s22*s31-s21*s32)*s43+(s21*s33-s23*s31)*s42+(s23*s32-s22*s33)*s41)+((s22*s31-s21*s32)*s44+(s21*s34-s24*s31)*s42+(s24*s32-s22*s34)*s41)*s53+((s21*s33-s23*s31)*s44+(s24*s31-s21*s34)*s43+(s23*s34-s24*s33)*s41)*s52+((s23*s32-s22*s33)*s44+(s22*s34-s24*s32)*s43+(s24*s33-s23*s34)*s42)*s51;
		
					w2 = m1*(((s32-s12)*s43+(s13-s33)*s42+s12*s33-s13*s32)*s54+((s12-s32)*s44+(s34-s14)*s42-s12*s34+s14*s32)*s53+((s33-s13)*s44+(s14-s34)*s43+s13*s34-s14*s33)*s52+(s13*s32-s12*s33)*s44+(s12*s34-s14*s32)*s43+(s14*s33-s13*s34)*s42)+m2*(((s11-s31)*s43+(s33-s13)*s41-s11*s33+s13*s31)*s54+((s31-s11)*s44+(s14-s34)*s41+s11*s34-s14*s31)*s53+((s13-s33)*s44+(s34-s14)*s43-s13*s34+s14*s33)*s51+(s11*s33-s13*s31)*s44+(s14*s31-s11*s34)*s43+(s13*s34-s14*s33)*s41)+m3*(((s31-s11)*s42+(s12-s32)*s41+s11*s32-s12*s31)*s54+((s11-s31)*s44+(s34-s14)*s41-s11*s34+s14*s31)*s52+((s32-s12)*s44+(s14-s34)*s42+s12*s34-s14*s32)*s51+(s12*s31-s11*s32)*s44+(s11*s34-s14*s31)*s42+(s14*s32-s12*s34)*s41)+((s12*s31-s11*s32)*s43+(s11*s33-s13*s31)*s42+(s13*s32-s12*s33)*s41)*s54+m4*(((s11-s31)*s42+(s32-s12)*s41-s11*s32+s12*s31)*s53+((s31-s11)*s43+(s13-s33)*s41+s11*s33-s13*s31)*s52+((s12-s32)*s43+(s33-s13)*s42-s12*s33+s13*s32)*s51+(s11*s32-s12*s31)*s43+(s13*s31-s11*s33)*s42+(s12*s33-s13*s32)*s41)+((s11*s32-s12*s31)*s44+(s14*s31-s11*s34)*s42+(s12*s34-s14*s32)*s41)*s53+((s13*s31-s11*s33)*s44+(s11*s34-s14*s31)*s43+(s14*s33-s13*s34)*s41)*s52+((s12*s33-s13*s32)*s44+(s14*s32-s12*s34)*s43+(s13*s34-s14*s33)*s42)*s51;

					w3 = m1*(((s12-s22)*s43+(s23-s13)*s42-s12*s23+s13*s22)*s54+((s22-s12)*s44+(s14-s24)*s42+s12*s24-s14*s22)*s53+((s13-s23)*s44+(s24-s14)*s43-s13*s24+s14*s23)*s52+(s12*s23-s13*s22)*s44+(s14*s22-s12*s24)*s43+(s13*s24-s14*s23)*s42)+m2*(((s21-s11)*s43+(s13-s23)*s41+s11*s23-s13*s21)*s54+((s11-s21)*s44+(s24-s14)*s41-s11*s24+s14*s21)*s53+((s23-s13)*s44+(s14-s24)*s43+s13*s24-s14*s23)*s51+(s13*s21-s11*s23)*s44+(s11*s24-s14*s21)*s43+(s14*s23-s13*s24)*s41)+m3*(((s11-s21)*s42+(s22-s12)*s41-s11*s22+s12*s21)*s54+((s21-s11)*s44+(s14-s24)*s41+s11*s24-s14*s21)*s52+((s12-s22)*s44+(s24-s14)*s42-s12*s24+s14*s22)*s51+(s11*s22-s12*s21)*s44+(s14*s21-s11*s24)*s42+(s12*s24-s14*s22)*s41)+((s11*s22-s12*s21)*s43+(s13*s21-s11*s23)*s42+(s12*s23-s13*s22)*s41)*s54+m4*(((s21-s11)*s42+(s12-s22)*s41+s11*s22-s12*s21)*s53+((s11-s21)*s43+(s23-s13)*s41-s11*s23+s13*s21)*s52+((s22-s12)*s43+(s13-s23)*s42+s12*s23-s13*s22)*s51+(s12*s21-s11*s22)*s43+(s11*s23-s13*s21)*s42+(s13*s22-s12*s23)*s41)+((s12*s21-s11*s22)*s44+(s11*s24-s14*s21)*s42+(s14*s22-s12*s24)*s41)*s53+((s11*s23-s13*s21)*s44+(s14*s21-s11*s24)*s43+(s13*s24-s14*s23)*s41)*s52+((s13*s22-s12*s23)*s44+(s12*s24-s14*s22)*s43+(s14*s23-s13*s24)*s42)*s51;
		
					w4 = m1*(((s22-s12)*s33+(s13-s23)*s32+s12*s23-s13*s22)*s54+((s12-s22)*s34+(s24-s14)*s32-s12*s24+s14*s22)*s53+((s23-s13)*s34+(s14-s24)*s33+s13*s24-s14*s23)*s52+(s13*s22-s12*s23)*s34+(s12*s24-s14*s22)*s33+(s14*s23-s13*s24)*s32)+m2*(((s11-s21)*s33+(s23-s13)*s31-s11*s23+s13*s21)*s54+((s21-s11)*s34+(s14-s24)*s31+s11*s24-s14*s21)*s53+((s13-s23)*s34+(s24-s14)*s33-s13*s24+s14*s23)*s51+(s11*s23-s13*s21)*s34+(s14*s21-s11*s24)*s33+(s13*s24-s14*s23)*s31)+m3*(((s21-s11)*s32+(s12-s22)*s31+s11*s22-s12*s21)*s54+((s11-s21)*s34+(s24-s14)*s31-s11*s24+s14*s21)*s52+((s22-s12)*s34+(s14-s24)*s32+s12*s24-s14*s22)*s51+(s12*s21-s11*s22)*s34+(s11*s24-s14*s21)*s32+(s14*s22-s12*s24)*s31)+((s12*s21-s11*s22)*s33+(s11*s23-s13*s21)*s32+(s13*s22-s12*s23)*s31)*s54+m4*(((s11-s21)*s32+(s22-s12)*s31-s11*s22+s12*s21)*s53+((s21-s11)*s33+(s13-s23)*s31+s11*s23-s13*s21)*s52+((s12-s22)*s33+(s23-s13)*s32-s12*s23+s13*s22)*s51+(s11*s22-s12*s21)*s33+(s13*s21-s11*s23)*s32+(s12*s23-s13*s22)*s31)+((s11*s22-s12*s21)*s34+(s14*s21-s11*s24)*s32+(s12*s24-s14*s22)*s31)*s53+((s13*s21-s11*s23)*s34+(s11*s24-s14*s21)*s33+(s14*s23-s13*s24)*s31)*s52+((s12*s23-s13*s22)*s34+(s14*s22-s12*s24)*s33+(s13*s24-s14*s23)*s32)*s51;
		
					w5 = m1*(((s12-s22)*s33+(s23-s13)*s32-s12*s23+s13*s22)*s44+((s22-s12)*s34+(s14-s24)*s32+s12*s24-s14*s22)*s43+((s13-s23)*s34+(s24-s14)*s33-s13*s24+s14*s23)*s42+(s12*s23-s13*s22)*s34+(s14*s22-s12*s24)*s33+(s13*s24-s14*s23)*s32)+m2*(((s21-s11)*s33+(s13-s23)*s31+s11*s23-s13*s21)*s44+((s11-s21)*s34+(s24-s14)*s31-s11*s24+s14*s21)*s43+((s23-s13)*s34+(s14-s24)*s33+s13*s24-s14*s23)*s41+(s13*s21-s11*s23)*s34+(s11*s24-s14*s21)*s33+(s14*s23-s13*s24)*s31)+m3*(((s11-s21)*s32+(s22-s12)*s31-s11*s22+s12*s21)*s44+((s21-s11)*s34+(s14-s24)*s31+s11*s24-s14*s21)*s42+((s12-s22)*s34+(s24-s14)*s32-s12*s24+s14*s22)*s41+(s11*s22-s12*s21)*s34+(s14*s21-s11*s24)*s32+(s12*s24-s14*s22)*s31)+((s11*s22-s12*s21)*s33+(s13*s21-s11*s23)*s32+(s12*s23-s13*s22)*s31)*s44+m4*(((s21-s11)*s32+(s12-s22)*s31+s11*s22-s12*s21)*s43+((s11-s21)*s33+(s23-s13)*s31-s11*s23+s13*s21)*s42+((s22-s12)*s33+(s13-s23)*s32+s12*s23-s13*s22)*s41+(s12*s21-s11*s22)*s33+(s11*s23-s13*s21)*s32+(s13*s22-s12*s23)*s31)+((s12*s21-s11*s22)*s34+(s11*s24-s14*s21)*s32+(s14*s22-s12*s24)*s31)*s43+((s11*s23-s13*s21)*s34+(s14*s21-s11*s24)*s33+(s13*s24-s14*s23)*s31)*s42+((s13*s22-s12*s23)*s34+(s12*s24-s14*s22)*s33+(s14*s23-s13*s24)*s32)*s41;

					best.w[0] = w1 / det;
					best.w[1] = w2 / det;
					best.w[2] = w3 / det;
					best.w[3] = w4 / det;
					best.w[4] = w5 / det;
				}
				else
				{
					best.w[0] = 0.0;
					best.w[1] = 0.0;
					best.w[2] = 0.0;
					best.w[3] = 0.0;
					best.w[4] = 0.0;
				}
			}
			
			// print average and standar deviation values in the top
			as<CharacterVector>(myList[0])[i*iter+ii] = i + 1;
			as<CharacterVector>(myList[1])[i*iter+ii] = 1.0;
			for( j = 0 ; j < nsource ; j++ )
			{
				as<CharacterVector>(myList[2+j])[i*iter+ii] = best.w[j];
			}
		}
	}
//	printf("\r                                            \r");

	Rcpp::DataFrame dfout(myList);
	return dfout;
}
