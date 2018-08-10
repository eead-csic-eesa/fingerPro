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
		return avg + gsl_ran_tdist(rng, (double)(n-1)) * dev / sqrt((double)n);
	}
	return avg;
}

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
	j = trials;
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
		
			if(ii == 0)
			{
				for( k = 0 ; k < vars ; k++ )
				{
					for( l = 0 ; l < nsource ; l++ )
					{
						source_c[l][k] = source[l][k];
					}
				}
			}
			// compute corrected sources
			else
			{
				for( k = 0 ; k < vars ; k++ )
				{
					for( l = 0 ; l < nsource ; l++ )
					{
						source_c[l][k] = correct(source[l][k], source_d[l][k], source_n[l]);
					}
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

				// measure error
				gof1 = 0.0;
				gof2 = 0.0;
				for( k = 0 ; k < vars ; k++ )
				{
					sum = 0.0;
					for( l = 0 ; l < nsource ; l++ )
					{
						// Non corrected mean
						//sum = sum + source[l][k] * try1[l];

						// Corrected mean
						sum = sum + source_c[l][k] * try1[l];
					}
					// GOF1
					gof1 = gof1 + fabs( point[i][k] - sum ) / ( max[k] - min[k] );

					// GOF2
					gof2 = gof2 + ( point[i][k] - sum ) * ( point[i][k] - sum ) / ( max[k] - min[k] ) / ( max[k] - min[k] );
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

			// print average and standar deviation values in the top
			as<CharacterVector>(myList[0])[i*iter+ii] = i + 1;
			as<CharacterVector>(myList[1])[i*iter+ii] = best.gof1;
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

