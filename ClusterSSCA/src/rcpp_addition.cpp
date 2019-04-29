#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
List RandomStart_cpp(const arma::vec& input_vector) {

  int i = sum(input_vector);
  int ncluster = input_vector.n_elem;
  arma::mat mem_cluster = arma::zeros(i, ncluster);
  arma::vec ind_cluster(i);

  int sign = 0;
  for (int j = 0; j < ncluster; j++) {
    for (int k = 0; k < input_vector(j); k++) {
      ind_cluster(sign + k) = j + 1;
    }
    sign = sign + input_vector(j);
  }
  arma::vec rd_ind_cluster = RcppArmadillo::sample(ind_cluster, i, FALSE);

  for (int l = 0; l < i; l++) {
    mem_cluster(l, (rd_ind_cluster(l) - 1)) = 1;
  }

  NumericVector tmp = wrap(rd_ind_cluster);
  tmp.attr("dim") = R_NilValue;

  List to_return = List::create(
    Named("mem") = mem_cluster,
    Named("ind") = tmp);
  return to_return;

}

//[[Rcpp::export]]
arma::mat MatrixCenter_cpp(const arma::mat& input_matrix, const int& center, const int& scale){
  /* the mean value of the colums are matrix with one row*/
  arma::mat variable_mean = arma::mean(input_matrix, 0);
  arma::mat variable_sd = arma::stddev(input_matrix, 0);
  int n_observation = input_matrix.n_rows;
  int n_variables = input_matrix.n_cols;
  arma::mat output_matrix(n_observation, n_variables);

  if (center) {
    for (int i = 0; i < n_variables; i ++){
      for (int j = 0; j < n_observation; j ++){
        output_matrix(j, i) = input_matrix(j, i) - variable_mean(0, i);
      }
    }
  }

  if (scale){
    for (int i = 0; i < n_variables; i ++){
      for (int j = 0; j < n_observation; j ++){
        output_matrix(j, i) = input_matrix(j, i) / variable_sd(0, i);
      }
    }
  }

  return (output_matrix);
}


//[[Rcpp::export]]
List sca_common_cpp(arma::mat data_con, int common) {

  /*
  known_sparse: whether the sparsity of the data is known (0 = unknown, otherwise the input is the sparsity of the data)
  data_con: concatenated data
  nvar: the vector indivates the number of variables in each block
  common: number of common components
  */

  /*
  int upper = 1e9;
  int all_var = data_con.n_cols;
  int all_mem = data_con.n_rows;
  */
  //Rprintf("entering sca\n");

  int all_var = data_con.n_cols;

  arma::mat data_con_t = MatrixCenter_cpp(data_con, 1, 0);
  arma::mat u;
  arma::vec s;
  arma::mat v;
  arma::svd_econ(u, s, v, data_con_t, "both", "dc");

  //Rprintf("middle sca1\n");

  arma::mat score_result = u.cols(0, (common - 1));
  arma::vec singular_result = s.subvec(0, (common - 1));
  arma::mat loading_raw_result = v.cols(0, (common - 1));

  arma::mat loading_result = arma::zeros(all_var, common);

  //Rprintf("middle sca2\n");
  for (int i = 0; i < common; i++) {
    for (int j = 0; j < all_var; j++) {
      loading_result(j, i) = loading_raw_result(j, i) * singular_result(i);
    }
  }
  arma::mat residule = data_con_t - score_result * trans(loading_result);
  arma::mat dore = residule % residule;
  double sum_residule = accu(dore);

  //Rprintf("exit sca\n");
  List to_return = List::create(
    Named("t") = score_result,
    Named("p") = loading_result,
    Named("l") = sum_residule);
  return to_return;
}

//[[Rcpp::export]]
List csca_cpp(arma::mat xc, arma::vec nvar, int nblock,
              int all_components, int ncluster, int iteration) {

  /*////////////////////////////////////////////////////////
  xc: concatenated data
  nvar: number of variables in each data block (the vector)
  nblock: number of blocks
  all_components: number of all components
  ncluster: number of clusters
  nrespondents: number of respondents
  iteration: number of maximum iteration
  ////////////////////////////////////////////////////////*/

  // set the upper bound of loss function
  //double upper = 100000;

  ///// set the minimum converge creterion
  ///// might be useful if we have converge scale
  // double converge_stage1 = 1e-2;
  //double converge_stage2 = 1e-5;

  // the aggregate level of information
  // int all_var = xc.n_cols;
  int all_mem = xc.n_rows;

  //set the upper bound of minimum loss
  double loss_min = 100000;
  int maxiter = 50;

  // collecte the total loss of each iteration
  arma::vec loss_all(iteration);

  // initialize the start
  arma::vec start(ncluster);
  for (int i = 0; i < ncluster - 1; i++) {
    start(i) = round(all_mem / ncluster);
  }
  start(ncluster - 1) = all_mem - (ncluster - 1) * round(all_mem / ncluster);

  // other memory that needs to be stated before the loop
  // arma::vec loss_iteration;
  // arma::mat cluster_yn(all_mem, ncluster);
  // arma::vec cluster_mem(all_mem);
  //////////////////////////////////////////////////////////////////////
  // step 2: perform the SSCA for each cluster and calculate the initial loss function
  ////////////////////////////////////////////////////////////////////////
  List current_data(ncluster);
  List current_loading(ncluster);
  List current_score(ncluster);
  arma::vec current_loss(ncluster);

  List temp_score(ncluster);
  List temp_loading(ncluster);
  arma::vec temp_loss(ncluster);
  List current_mem(ncluster);

  List opt_mem(ncluster);
  List opt_score(ncluster);
  List opt_loading(ncluster);
  arma::vec opt_cluster_mem;


  for (int v = 0; v < iteration; v++) {
    int flag = 0;
    double loss_tot = 0;
    //step 1: initialize the cluster assignment
    List random_start = RandomStart_cpp(start);

    // arma::mat cluster_yn = random_start[0];
    arma::vec cluster_mem = random_start[1];

    for (int i = 0; i < ncluster; i++) {

      arma::uvec pos = find(cluster_mem == (i + 1));
      arma::mat temp_data = xc.rows(pos);
      current_data[i] = temp_data;

      List ssca_results = sca_common_cpp(current_data[i], all_components);
      arma::mat c_l = ssca_results["p"];
      current_loading[i] = c_l;
      arma::mat c_s = ssca_results["t"];
      current_score[i] = c_s;

      double c_loss = as<double>(ssca_results[2]);
      current_loss(i) = c_loss;
      loss_tot = c_loss + loss_tot;
      current_mem[i] = pos;
    }

    // iteration
    // step 3: update cluster assignment for each observation
    int conv = 0;
    int iter = 0;
    int sign = 0;
    // note that loss_tot only update after one whole iteration, and loss_tot_ind update after every observation
    double loss_tot_ind = loss_tot;

    double loss_others = 0;
    double min_change = 100000;
    double loss_change = 100000;
    int mem_update = -1;

    while (conv == 0) {
      iter = iter + 1;

      // the sign of whethe there ar emembership exchanges
      sign = 0;
      // check the fit of a particular observation with a particular cluster
      for (int j = 0; j < all_mem; j++) {
        // set a impossible upper bound
        min_change = 100000;
        loss_change = 100000;
        mem_update = -1;
        //cluster_mem(j) = 1;
        for (int k = 0; k < ncluster; k++) {
          // check which cluster does the member belong to
          if (cluster_mem(j) != (k + 1)) {
            /////////the loss function with observation j
            //the data matrix with j
            arma::uvec current_m = current_mem[k];
            current_m.resize(current_m.size() + 1);
            current_m(current_m.size() - 1) = j;
            arma::mat new_data = xc.rows(current_m);

            // compute the new score matrix and loading matrix (use rational start, only compute one time)
            // compute the new loss values
            List new_ssca = sca_common_cpp(new_data, all_components);
            arma::mat temp_l = new_ssca[1];
            temp_loading[k] = temp_l;
            arma::mat temp_s = new_ssca[0];
            temp_score[k] = temp_s;
            double temp_los = as<double>(new_ssca[2]);
            temp_loss(k) = temp_los;

            loss_change = temp_loss(k) - current_loss(k);

            if (loss_change < min_change) {
              min_change = loss_change;
              // update the membership status of this certain observation
              mem_update = (k + 1);
              //Rprintf("check 12\n");
              //Rcout << min_change;
            }
          }

          if (cluster_mem(j) == (k + 1)) {
            //Rprintf("check 2\n");
            /////////// get the loss function without observation j
            //////  compute the data matrix without j
            // take the value from the list and convert it to the uvector
            arma::uvec current_m = current_mem[k];
            arma::uvec position_j = find(current_m != j);
            arma::uvec diff = current_m.elem(position_j);

            //arma::mat new_data = xc.rows(current_m);
            //arma::uvec current_m = current_mem[k];
            // using STD library to erase the member j in vector current_mem[k]
            //std::vector<int> current_m_std = arma::conv_to<std::vector<int>>::from(current_m);
             //Rprintf("check 3\n");
             //current_m_std.erase(current_m_std.begin() + j - 1);
            // convert back to uvector (which would be used in mat.rows)
            //arma::uvec diff = arma::conv_to<arma::uvec>::from(current_m_std);
            //Rprintf("check 4\n");
            arma::mat new_data_e = xc.rows(diff);

            // stop the loop when the new_data contains fewer than three rows
            if (new_data_e.n_rows < 3) {
              //set the flags to stop
              flag = 1;
              //Overall three loops need to be broken. Here we first break the for loop
              break;
            }
            //}

            ///// update both T and P
            // compute the new score matrix (rational start)
            List new_ssca = sca_common_cpp(new_data_e, all_components);
            arma::mat temp_l = new_ssca[1];
            temp_loading[k] = temp_l;
            arma::mat temp_s = new_ssca[0];
            temp_score[k] = temp_s;
            double temp_los = as<double>(new_ssca[2]);
            temp_loss(k) = temp_los;
            //Rcout << temp_loss(k);
            //Rcout << current_loss(k);

            loss_change = current_loss(k) - temp_loss(k);
            //Rcout << loss_change;
            //Rprintf("check 5\n");

            if (loss_change < min_change) {
              min_change = loss_change;
              // update the membership status of this certain observation
              mem_update = (k + 1);
              //Rprintf("check 12\n");
              //Rcout << min_change;
            }

            //Rprintf("where\n");
          }

          //Rprintf("are\n");
        }
        // break the loop j
        if (flag == 1) break;

        //update the cluster membership, loading matrix and the loss function (loss function only for indication) immediately
        //(only when the member is in differenet cluster)
        if (cluster_mem(j) != mem_update) {

          sign = 1;
          // update the cluster membership
          int c = cluster_mem(j) - 1;
          arma::uvec current_member = current_mem[c];
          arma::uvec position_c = find(current_member != j);
          arma::uvec diff_m = current_member.elem(position_c);
          //Rprintf("check 6\n");
          // using STD library to erase the member j in vector current_mem[k]
          //std::vector<int> current_member_std = arma::conv_to<std::vector<int>>::from(current_member);
          //current_member_std.erase(current_member_std.begin() + (j - 1));
          // convert back to uvector (which would be used in mat.rows)
          //arma::uvec diff_m = arma::conv_to<arma::uvec>::from(current_member_std);
          current_mem[c] = diff_m;
          //Rcout << diff_m;
          arma::uvec update_member = current_mem[(mem_update - 1)];
          update_member.resize(update_member.size() + 1);
          //Rprintf("check 8\n");
          update_member(update_member.size() - 1) = j;
          current_mem[(mem_update - 1)] = update_member;

          // loss function of other clusters
          loss_others = loss_tot_ind - current_loss(c) - current_loss((mem_update - 1));
          //Rprintf("check 9\n");

          ///////update the cluster which previously includes the current observation (exclude the current observation)
          //update the current loading
          arma::mat t_p = temp_loading[c];
          current_loading[c] = t_p;
          // update the current score
          arma::mat t_t = temp_score[c];
          current_loading[c] = t_t;
          // update the current loss
          current_loss(c) = temp_loss(c);
          //Rprintf("check 10\n");

          ////update the cluster which has the best fit with the current observation (include the current observation)
          // update the current loading
          arma::mat t_p_update = temp_loading[(mem_update - 1)];
          current_loading[(mem_update - 1)] = t_p_update;
          // updating the current score
          arma::mat t_t_update = temp_score[(mem_update - 1)];
          current_score[(mem_update - 1)] = t_t_update;
          // updating the current loss
          current_loss((mem_update - 1)) = temp_loss((mem_update - 1));
          //Rprintf("check 11\n");

          // update the overall membership vector
          cluster_mem(j) = mem_update;

          //update the loss function
          loss_tot_ind = loss_others + current_loss(c) + current_loss((mem_update - 1));
        }
      }
      if (flag == 1) break;

      // alternative way to see whether the cluster assignment has been updated.
      //double old_loss_tot = loss_tot;
      loss_tot = loss_tot_ind;
      //double loss_tot_change = old_loss_tot - loss_tot;

      if (sign == 0) {
        conv = 1;
      }
      if (iter == maxiter) {
        conv = 1;
      }
    }

    loss_all(v) = loss_tot;

    if (loss_tot < loss_min) {
      loss_min = loss_tot;
      opt_loading = current_loading;
      opt_score = current_score;
      opt_mem = current_mem;
      opt_cluster_mem = cluster_mem;
      //Rprintf("check 13\n");
    }
  }

  // convert the matrix loss_all to vector.
  NumericVector tmp = wrap(loss_all);
  tmp.attr("dim") = R_NilValue;

  //return the score matrices and loading matrices (not always necessary)

  /*
  arma::mat final_score(all_mem, all_components);
  List final_loading(ncluster);

  for (int t = 0; t < ncluster; t++) {
    arma::uvec opt_mem_ind = opt_mem[t];
    int length = opt_mem_ind.n_elem;
    for (int s = 0; s < length; s++) {
      arma::mat opt_score_ind = opt_score[t];
      arma::mat replace_mat = opt_score_ind.row(s);
      double tmp_value = opt_mem_ind(s);
      final_score.row(tmp_value) = replace_mat;
    }
    arma::mat opt_loading_ind = opt_loading[t];
    final_loading[t] = opt_loading_ind;
  }
  */

  //Rprintf("check 14\n");
  List to_return = List::create(
    Named("cluster_mem") = opt_cluster_mem,
    Named("loss") = loss_min,
    Named("loadings") = opt_loading,
    Named("scores") = opt_score
    );

  return to_return;
}

//[[Rcpp::export]]
List sparsedisco_cpp(arma::mat data_con, arma::vec nvar, int all_components, arma::mat p, arma::uvec fixed_zeros, int n_zeros) {

  /*
  data_con: the concatenated dataset
  nvar: a vector, each element indicates the number of variables in one data block
  all_components: total number of components (including common and distinctive components)
  p: the starting loading matrices
  fixed_zeros: numbers that indicate the positions of fixed zeros (i.e. structure-induced zeros)
  n_zeros: number of zeros in the loading matrices, which correspond to the assumed sparseness level

  */

  //number of variables and participants
  int j_tot = data_con.n_cols;
  int i_tot = data_con.n_rows;
  int p_row = p.n_rows;
  int p_col = p.n_cols;

  // loss functions
  double loss_p = 1000000;
  double update_loss;
  double stop = 1e-6;
  int maxiter = 100;

  int conv = 0;
  int iter = 0;

  // define all the variables which are used in the following algorithms
  arma::mat u;
  arma::vec s;
  arma::mat v;
  arma::mat t;
  arma::mat a_t;

  while (conv == 0) {
    iter = iter + 1;

    // conditional on P, re-estimate the score matrix T
    a_t = trans(data_con * p);

    arma::svd_econ(u, s, v, a_t, "both", "dc");
    t = v * trans(u);

    // update corresponds to each component (described in Gu & Van Deun 2016)
    // given the score matrix T, update the loading matrix P
    arma::mat pold = p;
    p = trans(data_con) * t;

    // impose zeros in loading matrices to represent the "distinctive components" (i.e. the variable "fixed_zeros" represent the position indexes)
    // Compare the absolute value of all non-zero elements in P, impose the smallest n_zeros elements to zero
    // first impose p to vectors (for the ease of later operations)
    arma::vec p_vec = arma::vectorise(p);
    // impose the structure-induced zeros (note that the positions of rcpp = (the positions of R) - 1)
    int fixed_zeros_n = fixed_zeros.size();
    arma::vec fixed_zeros_vec(fixed_zeros_n);
    fixed_zeros_vec.fill(0);
    p_vec.elem((fixed_zeros - 1)) = fixed_zeros_vec;

    // impose the sparse-induced zeros
    // sort the vector p_vec in ascending order
    arma::uvec sort_p_vec = sort_index(abs(p_vec));
    // select the smallest (fixed_zeros_n + n_zeros) elements (including zeros)
    arma::uvec first_p_vec = sort_p_vec.subvec(0, (fixed_zeros_n + n_zeros - 1));
    // converge these elements to zeros
    arma::vec sparse_zeros_vec(fixed_zeros_n + n_zeros);
    sparse_zeros_vec.fill(0);
    p_vec.elem(first_p_vec) = sparse_zeros_vec;

    // converge p back to the matrix
    //p.insert_cols(0, p_vec);
    p = reshape(p_vec, p_row, p_col);

    //update the current loss function
    arma::mat dev = data_con - t * trans(p);
    double ss = 0;
    for (int i = 0; i < i_tot; i++ ) {
      for (int j = 0; j < j_tot; j++ ) {
        ss = ss + (dev(i,j) * dev(i,j));
      }
    }

    update_loss = loss_p - ss;
    loss_p = ss;

  //stops
    if (maxiter == iter) {
      conv = 1;
    }
    if(update_loss < stop) {
      conv = 1;
    }
  }

    //return
    List to_return = List::create(
      Named("T") = t,
      Named("P") = p,
      Named("L") = loss_p,
      Named("iter") = iter);

    return to_return;
}

//[[Rcpp::export]]
List IntSparseSca_rational_full_cpp(arma::mat data_con, arma::vec nvar, int nblock, int all_components, arma::mat initial, arma::uvec fixed_zeros, double sparsity) {

  /*
  data_con: the concatenated dataset
  nvar: a vector, each element indicates the number of variables in one data block
  nblock: number of blocks in the full dataset
  all_components: total number of components (including common and distinctive components)
  fixed_zeros: numbers that indicate the positions of fixed zeros (i.e. structure-induced zeros)
  modes: the types of random generation mechanisms (1 = rational start (only rational start with SVD as starting decomposition);
  2 = semi-rational start (SVD starting partition with some additional noise); 3 = full (1+2))
  */

  int all_var = data_con.n_cols;
  int all_mem = data_con.n_rows;
  int fixed_zeros_n = fixed_zeros.size();
  int no_fixed_zeros = all_var * all_components - fixed_zeros_n;

  // centering the data
  arma::mat data_con_t = MatrixCenter_cpp(data_con, 1, 0);

  // total rational start
  arma::mat initial_p = initial;
  int zeros = round(no_fixed_zeros * sparsity);

  // find the solution (totally rational start)
  List new_ssca = sparsedisco_cpp(data_con_t, nvar, all_components, initial_p, fixed_zeros, zeros);
  arma::mat score_result = new_ssca[0];
  arma::mat loading_result = new_ssca[1];

  arma::mat residual = data_con_t - score_result * trans(loading_result);
  double ss = 0;
  for (int i = 0; i < all_mem; i++ ) {
    for (int j = 0; j < all_var; j++ ) {
      ss = ss + (residual(i,j) * residual(i,j));
    }
  }

  //return
  List to_return = List::create(
    Named("scores") = score_result,
    Named("loadings") = loading_result,
    Named("loss") = ss);

  return to_return;
}

//[[Rcpp::export]]
List IntSparseSca_random_cpp(arma::mat data_con, arma::vec nvar, int nblock, int all_components, arma::uvec fixed_zeros, double sparsity) {

  // similar to the
  double min_l = 1000000000;
  int all_var = data_con.n_cols;
  int all_mem = data_con.n_rows;
  int fixed_zeros_n = fixed_zeros.size();
  int no_fixed_zeros = all_var * all_components - fixed_zeros_n;

  // centering the data
  arma::mat data_con_t = MatrixCenter_cpp(data_con, 1, 0);
  int zeros = round(no_fixed_zeros * sparsity);

  // record the data
  List global(3);

 // the initial values of P: SVD + a number randomly simulated from N(0, 1) distribution
 // first compute the SVD base

  arma::mat u;
  arma::vec s;
  arma::mat v;
  arma::svd_econ(u, s, v, data_con_t, "both", "dc");


 arma::mat score_result = u.cols(0, (all_components - 1));
 arma::vec singular_result = s.subvec(0, (all_components - 1));
 arma::mat loading_raw_result = v.cols(0, (all_components - 1));

 arma::mat initial_p (all_var, all_components);

 for (int i = 0; i < all_components; i++) {
   for (int j = 0; j < all_var; j++) {
     initial_p(j, i) = loading_raw_result(j, i) * singular_result(i);
   }
 }

 // impose the structure-induced zeros (note that the positions of rcpp = (the positions of R) - 1)
 arma::vec fixed_zeros_vec(fixed_zeros_n);
 fixed_zeros_vec.fill(0);
 initial_p.elem((fixed_zeros - 1)) = fixed_zeros_vec;

 // the solution obtained from SVD base
 List new_ssca = sparsedisco_cpp(data_con_t, nvar, all_components, initial_p, fixed_zeros, zeros);
 double loss_result = as<double>(new_ssca[2]);
 if (loss_result < min_l){
   min_l = loss_result;
   global = new_ssca;
 }

 double solution_loss;
 arma::mat random_matrix_a;
 arma::mat start_p;

 // random starts (n = 150)
  for (int j = 0; j < 50; j++){
   random_matrix_a = 0 + 5 * arma::randn(all_var, all_components);
   start_p = initial_p + random_matrix_a;
   List solution =  sparsedisco_cpp(data_con_t, nvar, all_components, start_p, fixed_zeros, zeros);
   solution_loss = as<double>(solution[2]);
   if (solution_loss < min_l){
     min_l = solution_loss;
     global = solution;
   }
  }

  for (int j = 0; j < 50; j++){
    random_matrix_a = -5 + 5 * arma::randn(all_var, all_components);
    start_p = initial_p + random_matrix_a;
    List solution =  sparsedisco_cpp(data_con_t, nvar, all_components, start_p, fixed_zeros, zeros);
    solution_loss = as<double>(solution[2]);
    if (solution_loss < min_l){
      min_l = solution_loss;
      global = solution;
    }
  }

  for (int j = 0; j < 50; j++){
    random_matrix_a = 5 + 5 * arma::randn(all_var, all_components);
    start_p = initial_p + random_matrix_a;
    List solution =  sparsedisco_cpp(data_con_t, nvar, all_components, start_p, fixed_zeros, zeros);
    solution_loss = as<double>(solution[2]);
    if (solution_loss < min_l){
      min_l = solution_loss;
      global = solution;
    }
  }

  arma::mat score_result_prep = global[0];
  arma::mat loading_result_prep = global[1];
  score_result = score_result_prep;
  arma::mat loading_result = loading_result_prep;

  double ss = as<double>(global[2]);

  List to_return = List::create(
    Named("scores") = score_result,
    Named("loadings") = loading_result,
    Named("loss") = ss);

  return to_return;
}

//[[Rcpp::export]]
List cssca_quick_cpp(arma::mat data_con, arma::vec nvar, int nblock,
                   int common, arma::vec distinct, int ncluster, int nrespondents, double sparsity, arma::vec feed, double cutoff_prop) {
 //xc: concatenated data
 //nvar: number of variables in each data block
 //nblock: number of blocks
 //common: number of common components
 //dictinct: number of distinctive components for each data block
 //ncluster: number of clusters
 //memcluster: number of members in each cluster
 //sparsity: known sparsity
 //cutoff.prop: when the number of changing cluster memberships is smaller than this value,
 ////no random starting SSCA computation is needed

  // define several important variables
  double upper = 1000000000;
  int all_var = data_con.n_cols;
  int all_mem = data_con.n_rows;

  // number of all components (common + distinct)
  int n_distinct = 0;
  for (int i = 0; i < nblock; i++ ) {
    n_distinct = n_distinct + distinct(i);
  }
  int all_components = n_distinct + common;

  int flag = 0;
  // most iteration time for one single start
  int maxiter = 100;
  // loop index
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;

  // create the loading sparsity structure
  arma::vec distinct_index(n_distinct);
  arma::vec sum_var(nblock + 1);
  sum_var(0) = 0;
  int ini_var = 0;
  int count = 0;
  int i_distinct;
  for (i = 0; i < nblock; i++ ) {
     i_distinct = distinct(i); //1
     for (j = 0; j < i_distinct; j++){
       distinct_index(count) = (i + 1); // 1
       count = count + 1;// 1
     }
     ini_var = ini_var + nvar(i);
     sum_var(i + 1) = ini_var;
  }

  arma::uvec distinct_zeros;
  // initialize the distinct_zeros
  distinct_zeros = arma::linspace<arma::uvec>((all_var * common + sum_var(distinct_index(0) - 1) + 1), (all_var * common + sum_var(distinct_index(0))), (sum_var(distinct_index(0)) - sum_var(distinct_index(0) - 1)));
  for (i = 1; i < n_distinct; i++){
    arma::uvec temp_zeros = arma::linspace<arma::uvec>((all_var * (common + i) + sum_var(distinct_index(i) - 1) + 1), (all_var * (common + i) + sum_var(distinct_index(i))), (sum_var(distinct_index(i)) - sum_var(distinct_index(i) - 1)));
    distinct_zeros = join_cols(distinct_zeros, temp_zeros);
  }

  //set the upper bound of minimum loss
  double loss_min = upper;

  // the records
  arma::vec loss_iteration;
  double loss_tot = 0;

  //////////////////////////////////////////////////////////////////////
  // step 2: perform the SSCA for each cluster and calculate the initial loss function
  ////////////////////////////////////////////////////////////////////////


    //feed the cluster assignment
    arma::vec cluster_mem = feed;
    arma::vec cluster_mem_previous = feed;

    // record the current state
    List current_data(ncluster);
    List current_loading(ncluster);
    List current_score(ncluster);
    arma::vec current_loss(ncluster);
    List current_position(ncluster);
    List current_mem(ncluster);
    List temp_score(ncluster);
    List temp_loading(ncluster);
    arma::vec temp_loss(ncluster);

    // the additional records
    arma::vec current_mem_n(ncluster);
    arma::vec change_prop(ncluster);
    arma::vec change(ncluster);

    // the global results
    List opt_mem(ncluster);
    List opt_score(ncluster);
    List opt_loading(ncluster);
    arma::vec opt_cluster_mem;

    // the flags
    int conv = 0;
    int iter = 0;
    // the flags for checking whether any membershi has changed
    int sign = 0;

    // declare all the possible values
    double loss_change;

    while (conv == 0){
      iter = iter + 1;
      // the loss total differ in each iteration
      loss_tot = 0;
      // step 1: check the current number of memers in each cluster
      for (i = 0; i < ncluster; i++){
        // note the differences in index for C++ and R
        // only in the cluster_mem, we use the previous notation
        arma::uvec pos = find(cluster_mem == (i + 1));
        current_mem[i] = pos;
        current_mem_n(i) = pos.size();
        if (iter == 1){
          change(i) = current_mem_n(i);
        }
        change_prop(i) = change(i) / current_mem_n(i);

        // step 2: for each cluster, perform the initial full SSCA (it will be performed in each and every iteration)
        arma::mat temp_data = data_con.rows(pos);
        current_data[i] = temp_data;

        // some check (trial and fail)
        arma::vec test = arma::vectorise(temp_data);
        if (test.size() == 0){
          flag = 1;
          break;
        }
        if (arma::is_finite(test) == 0){
          flag = 1;
          break;
        }
        if (arma::is_finite(temp_data.n_rows) == 0 ){
          flag = 1;
          break;
        }
        if (temp_data.n_rows < (all_components + 1) ){
          flag = 1;
          break;
        }

        // special treatment for the first iteration
        if (iter == 1){
          List ssca_results = IntSparseSca_random_cpp(temp_data, nvar, nblock, all_components, distinct_zeros, sparsity);
          double loss_temp = as<double>(ssca_results[2]);
          loss_tot = loss_tot + loss_temp;
          current_loss(i) = loss_temp;
          arma::mat loading_temp = ssca_results[1];
          current_loading[i] = loading_temp;
          arma::mat score_temp = ssca_results[0];
          current_score[i] = score_temp;
          arma::vec loading_temp_vec = arma::vectorise(loading_temp);
          arma::uvec temp_pos = find(loading_temp_vec == 0);
          current_position[i] = temp_pos;
        }

        if (iter != 1){
          if (change_prop(i) > cutoff_prop){
            List ssca_results = IntSparseSca_random_cpp(temp_data, nvar, nblock, all_components, distinct_zeros, sparsity);
            double loss_temp = as<double>(ssca_results[2]);
            loss_tot = loss_tot + loss_temp;
            current_loss(i) = loss_temp;
            arma::mat loading_temp = ssca_results[1];
            current_loading[i] = loading_temp;
            arma::mat score_temp = ssca_results[0];
            current_score[i] = score_temp;
            arma::vec loading_temp_vec = arma::vectorise(loading_temp);
            arma::uvec temp_pos = find(loading_temp_vec == 0);
            current_position[i] = temp_pos;
          }
          else {
            arma::mat temp_loading = current_loading[i];
            List ssca_results = IntSparseSca_rational_full_cpp(temp_data, nvar, nblock, all_components, temp_loading, distinct_zeros, sparsity);
            double loss_temp = as<double>(ssca_results[2]);
            loss_tot = loss_tot + loss_temp;
            current_loss(i) = loss_temp;
            arma::mat loading_temp = ssca_results[1];
            current_loading[i] = loading_temp;
            arma::mat score_temp = ssca_results[0];
            current_score[i] = score_temp;
            arma::vec loading_temp_vec = arma::vectorise(loading_temp);
            arma::uvec temp_pos = find(loading_temp_vec == 0);
            current_position[i] = temp_pos;
          }
        }
      }

        if (flag == 1) {
          loss_tot = upper;
          break;
        }

        // records the number of changing observations and the future member assignments
        arma::vec change(ncluster);
        change.fill(0);

        arma::vec mem_update(all_mem);
        mem_update.fill(0);

        // add the loss_tot as the last element in loss_iteration
        loss_iteration.resize(loss_iteration.size() + 1);
        loss_iteration(loss_iteration.size() - 1) = loss_tot;
        if (iter > 1){
          if (loss_iteration(iter - 1) > loss_iteration(iter - 2)) {
            cluster_mem <- cluster_mem_previous;
            break;
          }
        }

         // Step 3: for each pair of the observation and the cluster, quickly update the corresponding loading matrices
         // with and without the observation. Update the cluster assignment of the observation so that the minimal loss functions could be
         //obtained

         for (i = 0; i < all_mem; i++){
           double min_change = upper;
           for (j = 0; j < ncluster; j++) {
             // when observation i belongs to cluster j (note the cluster index in cpp should be (cluster index in R - 1))
             if(cluster_mem(i) == (j + 1)){
               // get the data matrix without i
               // current_mem starting from 0
               arma::uvec current_m = current_mem[j];
               arma::uvec position_i = find(current_m != i);
               arma::uvec diff = current_m.elem(position_i);
               arma::mat  new_data_i = data_con.rows(diff);
               arma::mat  data_con_i = MatrixCenter_cpp(new_data_i, 1, 0);

               // update T and P for only one iteration
               arma::mat temp_loading = current_loading[j];
               arma::mat xp = trans(data_con_i * temp_loading);

               arma::mat u;
               arma::vec s;
               arma::mat v;
               arma::mat matrix_xi;

               arma::svd_econ(u, s, v, xp, "both", "dc");
               matrix_xi = v * trans(u);
               arma::mat p_new = trans(data_con_i) * matrix_xi;

               arma::uvec temp_position = current_position[j];


               // converge the matrix to the vector, impose the zeros on the vectors, then converge back to the matrix
               int temp_position_n = temp_position.size();
               arma::vec zero_vec(temp_position_n);
               zero_vec.fill(0);
               p_new.elem(temp_position) = zero_vec;

               arma::mat residual = data_con_i - matrix_xi * trans(p_new);

               int rows = data_con_i.n_rows;
               int cols = data_con_i.n_cols;

               double ss = 0;
               for (int k = 0; k < rows; k++ ) {
                 for (int l = 0; l < cols; l++ ) {
                   ss = ss + (residual(k,l) * residual(k,l));
                 }
               }

               double temp_loss = current_loss(j);
               loss_change = temp_loss - ss;
             }

             if(cluster_mem(i) != (j + 1)){

             //the data matrix without i
             arma::uvec current_m = current_mem[j];
             current_m.resize(current_m.size() + 1);
             current_m(current_m.size() - 1) = i;
             arma::mat new_data_i = data_con.rows(current_m);
             arma::mat  data_con_i = MatrixCenter_cpp(new_data_i, 1, 0);

             // update T and P for only one iteration
             arma::mat temp_loading = current_loading[j];
             arma::mat xp = trans(data_con_i * temp_loading);

             arma::mat u;
             arma::vec s;
             arma::mat v;
             arma::mat matrix_xi;

             arma::svd_econ(u, s, v, xp, "both", "dc");
             matrix_xi = v * trans(u);
             arma::mat p_new = trans(data_con_i) * matrix_xi;
             arma::uvec temp_position = current_position[j];

             // converge the matrix to the vector, impose the zeros on the vectors, then converge back to the matrix
             int temp_position_n = temp_position.size();
             arma::vec zero_vec(temp_position_n);
             zero_vec.fill(0);
             p_new.elem(temp_position) = zero_vec;

             arma::mat residual = data_con_i - matrix_xi * trans(p_new);

             int rows = data_con_i.n_rows;
             int cols = data_con_i.n_cols;

             double ss = 0;
             for (int k = 0; k < rows; k++ ) {
               for (int l = 0; l < cols; l++ ) {
                 ss = ss + (residual(k,l) * residual(k,l));
               }
             }

             double temp_loss = current_loss(j);
             loss_change = ss - temp_loss;
            }

            if (loss_change < min_change){
              min_change = loss_change;
              mem_update(i) = (j + 1);
            }
           }
         }
            //Step 4: formally update the cluster membership
            cluster_mem_previous = cluster_mem;

            for (i = 0; i < all_mem; i++){
                 if (cluster_mem(i) != mem_update(i)){
                   sign = 1;
                   int c = cluster_mem(i);
                   int c_new = mem_update(i);
                   cluster_mem(i) = mem_update(i);

                   // remove the observation
                   arma::uvec current_m = current_mem[(c - 1 )];
                   arma::uvec position_i = find(current_m != i);
                   arma::uvec diff = current_m.elem(position_i);
                   current_mem[(c - 1)] = diff;

                   // add the observation to the new cluster
                   arma::uvec current_m_1 = current_mem[(c_new - 1)];
                   current_m_1.resize(current_m_1.size() + 1);
                   current_m_1(current_m_1.size() - 1) = i;

                   change(c - 1) = change(c - 1) + 1;
                   change(c_new - 1) = change(c_new - 1) + 1;
                 }
            }

          if (sign == 0) {
            conv = 1;
          }
      }


    List to_return = List::create(
      Named("cluster_mem") = cluster_mem,
      Named("loss") = loss_tot,
      Named("scores") = current_score,
      Named("loadings") = current_loading
    );

    return to_return;
  }


