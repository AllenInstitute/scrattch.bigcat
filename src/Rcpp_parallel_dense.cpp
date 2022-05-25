//Rcpp_parallel_dense.cpp

// [[Rcpp::depends(beachmat)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include "beachmat3/beachmat.h"
#include <RcppParallel.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>
using namespace Rcpp;
using namespace RcppParallel;

///////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_means method1: iterate over only the non-zero elements in a column in a dense matrix
///////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_means_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// skip the column if it doesn't exist in clAll
		if (cl[col_i] < 1)
			continue;

		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		
		for (int row_i = 0; row_i < res.nrow(); row_i++)
		{
			res(row_i, mat_cluster_vector[col_i]) += *(indices + row_i);
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_means method2: iterate over only the non-zero elements in a column in a dense matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnMeanDense : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	const RVector<int> mat_cluster_vector;

	const RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source and destination
	ColumnMeanDense(beachmat::lin_matrix *matrix_ptr,
					 const IntegerVector mat_cluster_vector,
					 const IntegerVector cl,
					 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the row sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			if (cl[col_i] < 1)
				continue;

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < output.nrow(); row_i++)
			{
				output(row_i, mat_cluster_vector[col_i]) += *(indices + row_i);
			}
		}
	}
};

/*The function calls the ColumnMeanDense*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_means_RcppParallel_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean dense functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnMeanDense column_mean_dense(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_mean_dense, 20);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_present method1: iterate over only the non-zero elements in a column in a dense matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_present_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll, double lowth)
{

	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// skip the column if it doesn't exist in clAll
		if (cl[col_i] < 1)
			continue;

		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		

		for (int row_i = 0; row_i < res.nrow(); row_i++)
		{

			if (*(indices + row_i) < lowth)
				continue;

			res(row_i, mat_cluster_vector[col_i]) += 1.0;
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_present method2: iterate over only the non-zero elements in a column in a dense matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnPresentDense : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	RVector<int> mat_cluster_vector;

	RVector<int> cl;

	double lowth;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnPresentDense(beachmat::lin_matrix *matrix_ptr,
						IntegerVector mat_cluster_vector,
						IntegerVector cl,
						double lowth,
						NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), lowth(lowth), output(output) {}

	// get the column sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			if (cl[col_i] < 1)
				continue;

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			
			for (int row_i = 0; row_i < output.nrow(); row_i++)
			{

				if (*(indices + row_i) < lowth)
					continue;

				output(row_i, mat_cluster_vector[col_i]) += 1.0;
			}
		}
	}
};

/*The function calls the ColumnPresentSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_present_RcppParallel_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll, double lowth)
{

	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean dense functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnPresentDense column_present_dense(matrix_ptr, mat_cluster_vector, cl, lowth, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_present_dense, 20);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_means method1: iterate over only the non-zero elements in a column in a dense matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_means_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// skip the column if it doesn't exist in clAll
		if (cl[col_i] < 1)
			continue;

		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		

		for (int row_i = 0; row_i < res.nrow(); row_i++)
		{

			res(row_i, mat_cluster_vector[col_i]) += std::pow(*(indices + row_i),2.0);
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_means method2: iterate over only the non-zero elements in a column in a dense matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct ColumnSqrMeansDense : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	RVector<int> mat_cluster_vector;

	RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnSqrMeansDense(beachmat::lin_matrix *matrix_ptr,
						 IntegerVector mat_cluster_vector,
						 IntegerVector cl,
						 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the column sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			if (cl[col_i] < 1)
				continue;

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < output.nrow(); row_i++)
			{

				output(row_i, mat_cluster_vector[col_i]) += std::pow(*(indices + row_i), 2.0);
			}
		}
	}
};

/*The function calls the ColumnSqrMeansSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_means_RcppParallel_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean dense functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnSqrMeansDense column_sqr_means_dense(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_sqr_means_dense);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	for (int i = 0; i < res.ncol(); i++)
	{
		NumericMatrix::Column col = res(_, i);
		col = col / cluster_id_number[i];
	}

	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_median method1: iterate over only the non-zero elements in a column in a dense matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_medians_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name); //cl_pos is the position number of cl_all_cl_name in the count matrix

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash map: cluster_col_pos_map is the map for vector of position in each cluster*/
	std::unordered_map<int, std::vector<int>> cluster_col_pos_map;
	for (int i = 0; i < cl.size(); i++)
	{

		if (cl[i] < 1)
			continue;

		cluster_col_pos_map[cl[i]].push_back(i);
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster the column number in res matrix
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Loop through cluster_col_pos_map and create vector of vector*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), cluster_col_pos_map.size());
	std::vector<double> workspace(ptr->get_nrow());
	
	for (auto it = cluster_col_pos_map.begin(); it != cluster_col_pos_map.end(); ++it)
	{
		auto cluster_id = it->first;
		auto pos_vec = it->second;

		std::vector<std::vector<double>> sub_mat(ptr->get_nrow());
		for (auto pos_i : pos_vec)
		{
			// loop over col_i-th column of sub mat matrix
			auto indices = ptr->get_col(pos_i, workspace.data());
			

			for (int row_i = 0; row_i < sub_mat.size(); row_i++)
			{

				sub_mat[row_i].push_back(*(indices + row_i));
			}
		}

		size_t n = pos_vec.size() / 2;

		for (int i = 0; i < sub_mat.size(); i++)
		{
			sub_mat[i].resize(pos_vec.size(), 0);

			std::nth_element(sub_mat[i].begin(), sub_mat[i].begin() + n, sub_mat[i].end());
			double vn = sub_mat[i][n];

			if (pos_vec.size() % 2 == 1)
			{

				res(i, mat_cluster_vector[pos_vec[0]]) = vn;
			}
			else
			{

				std::nth_element(sub_mat[i].begin(), sub_mat[i].begin() + n - 1, sub_mat[i].end());
				res(i, mat_cluster_vector[pos_vec[0]]) = 0.5 * (vn + sub_mat[i][n - 1]);
			}
		}
	}

	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////
// get_cl_median method2: iterate over only the non-zero elements in a column in a dense matrix,using RcppParallel
//////////////////////////////////////////////////////////////////////////////////

struct ColumnMedianDense : public Worker
{
	// source matrix and maps
	std::vector<std::vector<double>> sub_mat;

	std::vector<int> mat_cluster_vector;

	std::vector<int> pos_vec;

	size_t n;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnMedianDense(std::vector<std::vector<double>> sub_mat,
					   std::vector<int> mat_cluster_vector,
					   std::vector<int> pos_vec,
					   size_t n,
					   NumericMatrix output)
		: sub_mat(sub_mat), mat_cluster_vector(mat_cluster_vector), pos_vec(pos_vec), n(n), output(output) {}

	// get the column sum of the given cols

	void operator()(std::size_t begin, std::size_t end)
	{

		for (int i = begin; i < end; i++)
		{
			sub_mat[i].resize(pos_vec.size(), 0);

			std::nth_element(sub_mat[i].begin(), sub_mat[i].begin() + n, sub_mat[i].end());
			double vn = sub_mat[i][n];

			if (pos_vec.size() % 2 == 1)
			{

				output(i, mat_cluster_vector[pos_vec[0]]) = vn;
			}
			else
			{

				std::nth_element(sub_mat[i].begin(), sub_mat[i].begin() + n - 1, sub_mat[i].end());
				output(i, mat_cluster_vector[pos_vec[0]]) = 0.5 * (vn + sub_mat[i][n - 1]);
			}
		}
	}
};

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_medians_RcppParallel_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name); //cl_pos is the position number of cl_all_cl_name in the count matrix

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash map: cluster_col_pos_map is the map for vector of position in each cluster*/
	std::unordered_map<int, std::vector<int>> cluster_col_pos_map;
	for (int i = 0; i < cl.size(); i++)
	{

		if (cl[i] < 1)
			continue;

		cluster_col_pos_map[cl[i]].push_back(i);
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster the column number in res matrix
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Loop through cluster_col_pos_map and create vector of vector*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), cluster_col_pos_map.size());
	std::vector<double> workspace(ptr->get_nrow());
	
	for (auto it = cluster_col_pos_map.begin(); it != cluster_col_pos_map.end(); ++it)
	{
		auto cluster_id = it->first;
		auto pos_vec = it->second;

		std::vector<std::vector<double>> sub_mat(ptr->get_nrow());
		for (auto pos_i : pos_vec)
		{
			// loop over col_i-th column of sub mat matrix
			auto indices = ptr->get_col(pos_i, workspace.data());
			

			for (int row_i = 0; row_i < sub_mat.size(); row_i++)
			{

				sub_mat[row_i].push_back(*(indices + row_i));
			}
		}

		size_t n = pos_vec.size() / 2;

		// column mean dense functor (pass input and output matrixes)

		ColumnMedianDense column_median_dense(sub_mat, mat_cluster_vector, pos_vec, n, res);

		// call parallelFor to do the work
		parallelFor(0, sub_mat.size(), column_median_dense);
	}

	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}





///////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sums method1: iterate over only the non-zero elements in a column in a dense matrix
///////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sums_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);


	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// skip the column if it doesn't exist in clAll
		if (cl[col_i] < 1)
			continue;

		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		

		for (int row_i = 0; row_i < res.nrow(); row_i++)
		{
			res(row_i, mat_cluster_vector[col_i]) += *(indices + row_i);
		}
	}

	
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sums method2: iterate over only the non-zero elements in a column in a dense matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnSumDense : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	const RVector<int> mat_cluster_vector;

	const RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source and destination
	ColumnSumDense(beachmat::lin_matrix *matrix_ptr,
					 const IntegerVector mat_cluster_vector,
					 const IntegerVector cl,
					 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the row sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			if (cl[col_i] < 1)
				continue;

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < output.nrow(); row_i++)
			{
				output(row_i, mat_cluster_vector[col_i]) += *(indices + row_i);
			}
		}
	}
};

/*The function calls the ColumnSumSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sums_RcppParallel_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{
	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector clAll_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(clAll_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(clAll_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(clAll_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	

	// column sum dense functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnSumDense column_sum_dense(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_sum_dense, 20);

	
	colnames(res) = colnames_vec;

	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}







//////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_sums method1: iterate over only the non-zero elements in a column in a dense matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_sums_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	std::vector<int> cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++)
	{
		// skip the column if it doesn't exist in clAll
		if (cl[col_i] < 1)
			continue;

		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace.data());
		

		for (int row_i = 0; row_i < res.nrow(); row_i++)
		{

			res(row_i, mat_cluster_vector[col_i]) += std::pow(*(indices + row_i),2.0);
		}
	}

	// loop over columns of res matrix, and divide by number of cells in each cluster
	//for (int i = 0; i < res.ncol(); i++)
	//{
	//	NumericMatrix::Column col = res(_, i);
	//	col = col / cluster_id_number[i];
	//}
	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_cl_sq_means method2: iterate over only the non-zero elements in a column in a dense matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct ColumnSqrSumsDense : public Worker
{
	// source matrix and maps
	beachmat::lin_matrix *matrix_ptr;

	RVector<int> mat_cluster_vector;

	RVector<int> cl;
	// destination matrix
	RMatrix<double> output;

	// initialize with source, cl_id and destination
	ColumnSqrSumsDense(beachmat::lin_matrix *matrix_ptr,
						 IntegerVector mat_cluster_vector,
						 IntegerVector cl,
						 NumericMatrix output)
		: matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), cl(cl), output(output) {}

	// get the column sum of the given cols
	void operator()(std::size_t begin, std::size_t end)
	{
		std::vector<double> workspace(matrix_ptr->get_nrow());
		
		for (int col_i = begin; col_i < end; col_i++)
		{

			if (cl[col_i] < 1)
				continue;

			// loop over col_i-th column of mat matrix
			auto indices = matrix_ptr->get_col(col_i, workspace.data());
			

			for (int row_i = 0; row_i < output.nrow(); row_i++)
			{

				output(row_i, mat_cluster_vector[col_i]) += std::pow(*(indices + row_i), 2.0);
			}
		}
	}
};

/*The function calls the ColumnSqrSumsSparse*/
//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_cl_sqr_sums_RcppParallel_dense(Rcpp::RObject mat, Rcpp::IntegerVector clAll)
{

	// call R function to find positions of each cell in matrix
	Function colnames_f("colnames");
	CharacterVector mat_name = colnames_f(mat);

	Function names_f("names");
	CharacterVector cl_all_cl_name = names_f(clAll);

	Function match_f("match");
	IntegerVector cl_pos = match_f(cl_all_cl_name, mat_name);

	CharacterVector levs_all = clAll.attr("levels");
	std::unordered_map<std::string, int> cluster_col_name_id_map;
	std::unordered_map<std::string, std::string> cluster_col_name_lev_map;

	for (int i = 0; i < clAll.length(); i++)
	{
		cluster_col_name_id_map[std::string(cl_all_cl_name[i])] = clAll[i];
		cluster_col_name_lev_map[std::string(cl_all_cl_name[i])] = std::string(levs_all[clAll[i] - 1]);
	}
	IntegerVector cl(mat_name.length(), 0);
	CharacterVector cl_colnames_vec(mat_name.length());
	for (int i = 0; i < cl_pos.length(); i++)
	{
		cl[cl_pos[i] - 1] = cluster_col_name_id_map[std::string(mat_name[cl_pos[i] - 1])];
		cl_colnames_vec[cl_pos[i] - 1] = cluster_col_name_lev_map[std::string(mat_name[cl_pos[i] - 1])];
	}

	/* Set up hash maps for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	   cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	IntegerVector mat_cluster_vector(cl.size(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector colnames_vec;
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end())
		{
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(cl_colnames_vec[i]);
		}
		else
		{
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}

	/* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for (int i = 0; i < cl.size(); i++)
	{
		if (cl[i] < 1)
			continue;
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}

	/* Loop though dense matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the dense matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
	std::vector<double> workspace(ptr->get_nrow());
	
	// column mean dense functor (pass input and output matrixes)
	auto matrix_ptr = ptr.get();
	ColumnSqrSumsDense column_sqr_sums_dense(matrix_ptr, mat_cluster_vector, cl, res);

	// call parallelFor to do the work
	parallelFor(0, ptr->get_ncol(), column_sqr_sums_dense);

	// loop over columns of res matrix, and divide by number of cells in each cluster
	//for (int i = 0; i < res.ncol(); i++)
	//{
	//	NumericMatrix::Column col = res(_, i);
	//	col = col / cluster_id_number[i];
	//}

	colnames(res) = colnames_vec;
	// call R function to read and assign rownames
	Function rownames_f("rownames");
	rownames(res) = rownames_f(mat);
	return res;
}
