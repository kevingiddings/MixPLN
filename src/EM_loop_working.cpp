#include <RcppArmadillo.h>
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Helper functions

// Function to take a std::vector as first argument, scalar as second, and add scalar to every entry
std::vector<long double> vec_add_const(const std::vector<long double> input_vec, const long double scalar)
{
    std::size_t num_elements = input_vec.size();

    std::vector<long double> return_vec(num_elements);

    for (int i = 0; i < num_elements; i++)
    {
        return_vec[i] = input_vec[i] + scalar;
    }
    return return_vec;
}

// Function to take an arma mat, a row number, and return that row in the matrix it as a std::vector
std::vector<long double> retrieve_row(const arma::mat input_mat, const int row_num)
{
    int num_cols = input_mat.n_cols;

    std::vector<long double> return_vec(num_cols);

    for (int i = 0; i < num_cols; i++)
    {
        long double val = input_mat(row_num, i);
        return_vec[i] = val;
    }
    return return_vec;
}

// Function to retireve max element of a vector
long double max_entry(const std::vector<long double> input_vec)
{
    std::size_t num_elements = input_vec.size();

    long double maximum = input_vec[0];

    for (int i = 1; i < num_elements; i++)
    {
        if (input_vec[i] > maximum)
        {
            maximum = input_vec[i];
        }
    }
    return maximum;
}

// Function to add values of two vectors of same length
std::vector<long double> add_vectors(const std::vector<long double> input_vec_1, const std::vector<long double> input_vec_2)
{
    std::size_t num_elements_1 = input_vec_1.size();
    // std::size_t num_elements_2 = input_vec_2.size();

    std::vector<long double> return_vec(num_elements_1);

    for (int i = 0; i < num_elements_1; i++)
    {
        return_vec[i] = input_vec_1[i] + input_vec_2[i];
    }
    return return_vec;
}

// Function to vectorize exp function
std::vector<long double> vec_exp(const std::vector<long double> input_vec)
{
    std::size_t num_elements = input_vec.size();

    std::vector<long double> return_vec(num_elements);

    for (int i = 0; i < num_elements; i++)
    {
        long double val = exp(input_vec[i]);
        return_vec[i] = val;
    }

    return return_vec;
}

// Function to compute sum of a vector
long double vec_sum(const std::vector<long double> input_vec)
{
    // Determine the number of elements in the input vector
    std::size_t num_elements = input_vec.size();

    // Initialize the sum
    long double sum = 0;

    // Compute the sum
    for (std::size_t i = 0; i < num_elements; i++)
    {
        sum += input_vec[i];
    }

    return sum;
}

// Function to compute log of entries in vector
std::vector<long double> nat_log(const arma::vec input_vec)
{
    // Determine the number of elements in the input vector
    std::size_t num_elements = input_vec.size();

    // Create a vector to store the natural logarithms of the input vector elements
    std::vector<long double> result(num_elements);

    // Compute natural logarithms
    for (std::size_t i = 0; i < num_elements; i++)
    {
        long double log_i = log(input_vec(i));
        result[i] = log_i;
    }

    return result;
}

// Function to compute row sums of matrix
std::vector<long double> row_sums(const arma::mat input_matrix)
{
    // Determine the number of rows and columns in the input matrix
    std::size_t num_rows = input_matrix.n_rows;
    std::size_t num_cols = input_matrix.n_cols;

    // Create a matrix to store the row sums
    std::vector<long double> result(num_rows);

    // Compute row sums
    for (int i = 0; i < num_rows; i++)
    {
        long double sum = 0;
        for (int j = 0; j < num_cols; j++)
        {
            sum += input_matrix(i, j);
        }
        result[i] = sum;
    }

    return result;
}

arma::mat divide_by_row_sums(const arma::mat input_mat)
{
    std::size_t num_rows = input_mat.n_rows;
    std::size_t num_cols = input_mat.n_cols;

    arma::mat result(num_rows, num_cols);
    std::vector<long double> row_sums_vec = row_sums(input_mat);

    // Divide each element of F by the corresponding row sum
    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            result(i, j) = input_mat(i, j) / row_sums_vec[i];
        }
    }

    return result;
}

arma::mat for_mu_sig(int i, int g, Rcpp::List iOsigO_g, arma::mat z, Rcpp::List O_list, arma::mat m_g_i)
{
    arma::mat O_list_i = Rcpp::as<mat>(O_list[i]);
    arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);

    arma::mat temp = z(i, g) * trans(O_list_i) * iOsigO_g_i;
    arma::vec num = temp * m_g_i;
    arma::mat den = temp * O_list_i;

    arma::mat result = join_rows(den, num);

    return result;
}

arma::mat for_sig(int i, int g, Rcpp::List iOsigO_g, arma::mat z, Rcpp::List O_list, arma::vec m_g_i, Rcpp::List mu, arma::mat S_g_i)
{
    arma::mat O_list_i = Rcpp::as<mat>(O_list[i]);
    arma::vec mu_g = Rcpp::as<vec>(mu[g]);
    arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);

    arma::mat omega = (m_g_i - O_list_i * mu_g) * (m_g_i - O_list_i * mu_g).t() + S_g_i;
    arma::mat forsig1 = z(i, g) * trans(O_list_i) * iOsigO_g_i * omega * iOsigO_g_i * O_list_i;
    arma::mat forsig2 = z(i, g) * trans(O_list_i) * iOsigO_g_i * O_list_i;
    return -(forsig1 - forsig2);
}

arma::vec nan_omit(arma::vec x)
{
    arma::vec y = x;
    for (int i = 0; i < x.n_elem; i++)
    {
        if (std::isnan(x[i]))
        {
            y[i] = 0.0L;
        }
    }
    return y;
}

arma::vec lfactorial(arma::vec x)
{
    arma::vec y = x;
    for (int i = 0; i < x.n_elem; i++)
    {
        y[i] = lgamma(x[i] + 1.0L);
    }
    return y;
}

long double abs_ld(const long double input){
    if(input < 0){
        return -input;
    }else{
        return input;
    }
}

// [[Rcpp::export]]
Rcpp::List EM_loop(Rcpp::List &start, arma::vec &lib_mat, Rcpp::List &m, Rcpp::List mu, Rcpp::List sigma,
                   Rcpp::List isigma, Rcpp::List iOsigO, Rcpp::List &S, Rcpp::List &O_list, arma::mat &Y,
                   arma::mat z, arma::mat &O_mat, Rcpp::List pi_g, const int &d, const int &N, const int &G)
{
    // std::ofstream myfile;
    // myfile.open("Dang_new_output9.txt");
    // myfile << std::fixed << std::setprecision(10);
    Rcpp::List sigma_new(G);
    Rcpp::List GX(G);
    Rcpp::List dGX(G);
    Rcpp::List z_S(G);

    double step = 0.001L;
    double checks = 0L;
    int it = 1L;
    int it_max = 1000L;
    std::vector<long double> aloglik = {0.0L, 0.0L, 0.0L};
    aloglik.resize(N);
    std::vector<long double> loglik;


    // NOV 14 update: moved z(N,G) to here, possibly was being reset every iteration?
    arma::mat z_loop(N,G);
    z_loop = z;  // Ensure z is initialized with the correct dimensions

    std::string exit_code = "UNKNOWN STOPPING; POSSIBLE BUG";

    Rcpp::Rcout << "Starting EM loop" << std::endl;
    while (checks == 0)
    {
        for (int g = 0; g < G; g++)
        {
            // Initialize or update GX, dGX, z_S, etc.
            // Rcpp::List GX_g(N), dGX_g(N), z_S_g(N), m_g(N), S_g(N), iOsigO_g(N), start_g(N);
            Rcpp::List GX_g(N);
            Rcpp::List dGX_g(N);
            Rcpp::List z_S_g(N);
            Rcpp::List m_g = Rcpp::as<Rcpp::List>(m[g]);
            Rcpp::List S_g = Rcpp::as<Rcpp::List>(S[g]);
            Rcpp::List start_g = Rcpp::as<Rcpp::List>(start[g]);
            Rcpp::List iOsigO_g = Rcpp::as<Rcpp::List>(iOsigO[g]);

            arma::vec mu_g = Rcpp::as<vec>(mu[g]);
            arma::mat sigma_g = Rcpp::as<mat>(sigma[g]);
            arma::mat isigma_g = Rcpp::as<mat>(isigma[g]);

            for (int i = 0; i < N; i++)
            {
                arma::mat O_list_i = Rcpp::as<mat>(O_list[i]);
                arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);
                arma::vec start_g_i = Rcpp::as<vec>(start_g[i]);
                arma::mat S_g_i = Rcpp::as<mat>(S_g[i]);

                arma::mat dGX_g_i = diagmat(arma::exp(O_list_i * vectorise(arma::log(lib_mat)) + start_g_i) + 0.5 * S_g_i.diag()) + iOsigO_g_i;
                S_g_i = inv(dGX_g_i);
                arma::mat z_S_g_i = z_loop(i, g) * S_g_i;
                arma::vec GX_g_i = O_list_i * vectorise(Y.row(i)) -
                                   vectorise(arma::exp(start_g_i + O_list_i * vectorise(arma::log(lib_mat)) + 0.5 * S_g_i.diag())) -
                                   iOsigO_g_i * (start_g_i - O_list_i * vectorise(mu_g));
                arma::vec m_g_i = start_g_i + S_g_i * vectorise(GX_g_i);
                // save varaibles to the list
                GX_g[i] = GX_g_i;
                dGX_g[i] = dGX_g_i;
                z_S_g[i] = z_S_g_i;
                m_g[i] = m_g_i;
                S_g[i] = S_g_i;
            }
            start[g] = m[g];

            // Update mu and sigma
            arma::mat temp3 = zeros(d, d + 1);
            // loop over all values of i, 1 to N, and add results of for_mu_sig to temp3
            for (int i = 0; i < N; i++)
            {
                arma::mat for_mu_sig_i = for_mu_sig(i, g, iOsigO_g, z_loop, O_list, m_g[i]);
                temp3 += for_mu_sig_i;
            }

            mu_g = inv(temp3.submat(0, 0, d - 1, d - 1)) * temp3.col(d);
            arma::mat gr = zeros(d, d);
            for (int i = 0; i < N; i++)
            {
                arma::mat for_sig_i = for_sig(i, g, iOsigO_g, z_loop, O_list, m_g[i], mu, S_g[i]);
                gr += for_sig_i;
            }

            arma::mat sigma_new_g = sigma_g - step * gr;
            sigma_new[g] = sigma_new_g;
            
            
            // Save varaibles to list
            //Rcpp::Rcout << "Saving variables..." << std::endl;
            GX[g] = GX_g;
            dGX[g] = dGX_g;
            z_S[g] = z_S_g;
            m[g] = m_g;
            S[g] = S_g;
            start[g] = start_g;
            iOsigO[g] = iOsigO_g;
        }


        //Boolean vector to indicate which new sigma matricies are semi positive definite
        std::vector<bool> Sigmas_PD(G);
        for(int g = 0; g < G; g++){
            arma::mat sigma_new_g = sigma_new[g];
            Sigmas_PD[g] = sigma_new_g.is_sympd();
        }
        
        //Check if all matricies are SPD
        bool all_SPD;
        for(int g = 0; g < G; g++){
            if(Sigmas_PD[g]){
                if(g == G - 1){
                    all_SPD = false;
                }
                continue;;
            }else{
                all_SPD = false;
                break;
            }
        }

        //if all new sigmas are SPD, then we update all the sigmas to the new values
        for(int g = 0; g < G; g++){
            if(all_SPD){
                sigma[g] = sigma_new[g];
            }
            arma::mat sigma_g = sigma[g];
            isigma[g] = inv(sigma_g);
        }

        pi_g = sum(z_loop, 0) / N;
        arma::mat lib_mat_full(N, lib_mat.size());
        // Fill each row of the result matrix with lib_mat
        for (int i = 0; i < N; i++)
        {
            lib_mat_full.row(i) = lib_mat.t(); // Transpose to convert vector to row vector
        }
        //Rcpp::Rcout << "Calculating iOsigO..." << std::endl;
        for (int g = 0; g < G; g++)
        {
            Rcpp::List iOsigO_g(N);
            iOsigO_g = iOsigO[g];
            arma::mat sigma_g = Rcpp::as<mat>(sigma[g]);
            for (int i = 0; i < N; i++)
            {
                arma::mat O_list_i = Rcpp::as<mat>(O_list[i]);
                iOsigO_g[i] = inv(O_list_i * sigma_g * trans(O_list_i));
            }
            iOsigO[g] = iOsigO_g;
        }

        arma::mat F(N, G);
        arma::mat F_raw(N,G);
        //Rcpp::Rcout << "Calculating F(i,g) matrix..." << std::endl;
        for (int g = 0; g < G; g++)
        {
            Rcpp::List iOsigO_g = Rcpp::as<Rcpp::List>(iOsigO[g]);
            Rcpp::List S_g = Rcpp::as<Rcpp::List>(S[g]);
            Rcpp::List m_g = Rcpp::as<Rcpp::List>(m[g]);
            for (int i = 0; i < N; i++)
            {
                arma::mat S_g_i = Rcpp::as<mat>(S_g[i]);
                arma::mat O_list_i = Rcpp::as<mat>(O_list[i]);
                arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);
                arma::vec m_g_i = Rcpp::as<vec>(m_g[i]);
                arma::vec mu_g = Rcpp::as<vec>(mu[g]);

                long double t_1_1 = log(det(S_g_i));
                long double term_1 = 0.5L * t_1_1;

                arma::rowvec t_2_1 = conv_to<rowvec>::from(m_g_i - O_list_i * vectorise(mu_g));
                arma::vec t_2_3 = m_g_i - O_list_i * vectorise(mu_g);
                long double t_2_2 = as_scalar(t_2_1 * iOsigO_g_i * vectorise(t_2_3));
                long double term_2 = -0.5L * t_2_2;

                long double t_3_1 = arma::trace(iOsigO_g_i * S_g_i);
                long double term_3 = -t_3_1;

                long double t_4_1 = log(det(iOsigO_g_i));
                long double term_4 = 0.5L * t_4_1;

                long double t_5_1 = arma::sum(nan_omit(vectorise(O_mat.row(i))));
                long double term_5 = 0.5L * t_5_1;

                arma::rowvec m_g_i_row = conv_to<rowvec>::from(m_g_i);
                long double term_6 = as_scalar(m_g_i_row * O_list_i * vectorise(Y.row(i)));

                long double t_7_1 = arma::sum(vectorise(arma::exp(m_g_i + 0.5 * arma::diagvec(S_g_i))) + lfactorial(O_list_i * vectorise(Y.row(i))));
                long double term_7 = -t_7_1;

                long double total = term_1 + term_2 + term_3 + term_4 + term_5 + term_6 + term_7;

                // Looks silly, but this is required due to limitations of double precision in arma objects
                // Need long double precision due to F_i_g being very close to 0.
                long double e_total = expl(total);
                long double pi_g_g = pi_g[g];
                long double F_i_g = pi_g_g * e_total;
                F(i, g) = F_i_g;
                F_raw(i, g) = total;

                // myfile << "Iteration: " << it << std::endl;
                // myfile << "Group: " << g + 1<< std::endl;
                // myfile << "Row: " << i + 1<< std::endl;
                // myfile << "Term 1: " << term_1 << std::endl;
                // myfile << "Term 2: " << term_2 << std::endl;
                // myfile << "Term 3: " << term_3 << std::endl;
                // myfile << "Term 4: " << term_4 << std::endl;
                // myfile << "Term 5: " << term_5 << std::endl;
                // myfile << "Term 6: " << term_6 << std::endl;
                // myfile << "Term 7: " << term_7 << std::endl;
                // myfile << "Total: " << total << std::endl;
                // myfile << "F_i_g: " << F_i_g << std::endl;
            }
        }

        loglik.resize(it);
        // new loglik calculation
        //Rcpp::Rcout << "Calculating Log Likelihood..." << std::endl;
        loglik[it - 1] = 0;
        for (int i = 0; i < N; i++) {
            try {
                //Rcpp::Rcout << "Iteration i = " << i << std::endl;

                // Original code with diagnostics
                std::vector<long double> x_i(G);
                //Rcpp::Rcout << "Initialized x_i vector for i = " << i << " of size G = " << G << std::endl;

                for (int g = 0; g < G; g++) {
                    //Rcpp::Rcout << "Processing g = " << g << std::endl;
                    long double pi_g_g = pi_g[g];  // Ensure pi_g[g] is within bounds
                    //Rcpp::Rcout << "pi_g_g = " << pi_g_g << std::endl;
                    x_i[g] = log(pi_g_g) + F_raw(i, g);
                }

                //Rcpp::Rcout << "x_i populated: " << x_i[0] << std::endl;  // Optionally print more elements of x_i

                // Calculate C (maximum entry)
                long double C = floor(max_entry(x_i));
                //Rcpp::Rcout << "max_entry C = " << C << std::endl;


                // Step 1: Apply vec_add_const(x_i, -C)
                std::vector<long double> x_i_adjusted = vec_add_const(x_i, -C);
                //Rcpp::Rcout << "x_i_adjusted: ";
                //for (auto val : x_i_adjusted) Rcpp::Rcout << val << " ";
                //Rcpp::Rcout << std::endl;

                // Step 2: Apply vec_exp to x_i_adjusted
                std::vector<long double> x_i_exp = vec_exp(x_i_adjusted);
                //Rcpp::Rcout << "x_i_exp: ";
                //for (auto val : x_i_exp) Rcpp::Rcout << val << " ";
                //Rcpp::Rcout << std::endl;

                // Step 3: Sum the elements of x_i_exp using vec_sum
                long double exp_sum = vec_sum(x_i_exp);
                //Rcpp::Rcout << "exp_sum = " << exp_sum << std::endl;

                // Step 4: Calculate the log of exp_sum
                long double log_exp_sum = log(exp_sum);
                //Rcpp::Rcout << "log(exp_sum) = " << log_exp_sum << std::endl;

                // Step 5: Calculate 'change' as C + log(exp_sum)
                long double change = C + log_exp_sum;
                //Rcpp::Rcout << "change = " << change << std::endl;

                loglik[it - 1] += change;
                //Rcpp::Rcout << "Updated loglik[" << it - 1 << "] = " << loglik[it - 1] << std::endl;

            } catch (std::exception &ex) {
                Rcpp::Rcout << "Error: log-like update failed at i = " << i << " with message: " << ex.what() << std::endl;
                return Rcpp::List();  // Stop execution in case of an error
            } catch (...) {
                Rcpp::Rcout << "Error: log-like update failed (unknown issue) at i = " << i << std::endl;
                return Rcpp::List();  // Catch any unforeseen exceptions
            }
        }

        //Rcpp::Rcout << "Log-likelihood calculation completed successfully." << std::endl;

        // new z_ig calculation
        //Rcpp::Rcout << "Calculating Z(i,g) matrix..." << std::endl;
        for (int i = 0; i < N; i++) {
            try {
                //Rcpp::Rcout << "Processing row i = " << i << std::endl;

                // Step 1: Initialize x_i and populate it with log(pi_g[g]) + F(i, g)
                std::vector<long double> x_i(G);
                for (int g = 0; g < G; g++) {
                    long double pi_g_g = pi_g[g];  // Access pi_g[g]
                    x_i[g] = log(pi_g_g) + F_raw(i, g);  // Compute log(pi_g[g]) + F(i, g)
                }
                //Rcpp::Rcout << "x_i: ";
                //for (auto val : x_i) Rcpp::Rcout << val << " ";
                //Rcpp::Rcout << std::endl;

                // Step 2: Compute the maximum entry C from x_i
                long double C = floor(max_entry(x_i));
                //Rcpp::Rcout << "max_entry C = " << C << std::endl;

                // Step 3: Adjust x_i by subtracting C (x_i_adjusted = x_i - C)
                std::vector<long double> x_i_adjusted = vec_add_const(x_i, -C);
                //Rcpp::Rcout << "x_i_adjusted: ";
                //for (auto val : x_i_adjusted) Rcpp::Rcout << val << " ";
                //Rcpp::Rcout << std::endl;

                // Step 4: Apply exp to each element of x_i_adjusted
                std::vector<long double> x_i_exp = vec_exp(x_i_adjusted);
                //Rcpp::Rcout << "x_i_exp: ";
                //for (auto val : x_i_exp) Rcpp::Rcout << val << " ";
                //Rcpp::Rcout << std::endl;

                // Step 5: Compute the sum of the exponentials
                long double exp_sum = vec_sum(x_i_exp);
                //Rcpp::Rcout << "exp_sum = " << exp_sum << std::endl;

                // Step 6: Compute the Z(i, g) values and assign them to the matrix z
                for (int g = 0; g < G; g++) {
                    long double z_ig = exp(x_i[g] - C - log(exp_sum));  // Calculate Z(i, g)
                    z_loop(i, g) = z_ig;  // Assign Z(i, g) to the matrix
                    //Rcpp::Rcout << "z(" << i << ", " << g << ") = " << z_ig << std::endl;
                }

            } catch (std::exception &ex) {
                Rcpp::Rcout << "Error: z_loop(i,g) update failed at i = " << i << ". Message: " 
                            << ex.what() << std::endl;
                return Rcpp::List();  // Stop execution in case of an error
            } catch (...) {
                Rcpp::Rcout << "Error: z_loop(i,g) update failed (unknown issue) at i = " << i << std::endl;
                return Rcpp::List();  // Catch any unforeseen exceptions
            }
        }

        // old z_ig calculation
        // arma::mat z = divide_by_row_sums(F);
        //Rcpp::Rcout << "Calculating Aitken acceleration..." << std::endl;
        if (it > 3)
        {
            if (loglik[it - 2] - loglik[it - 3] == 0)
            {
                exit_code = "Log Likelihood equal for two iterations";
                checks = 1;
            }
            else
            {
                long double l1 = loglik[it - 3];
                long double l2 = loglik[it - 2];
                long double l3 = loglik[it - 1];

                long double a;
                if (it < 5 || l2 > l1)
                {
                    a = (l3 - l2) / (l2 - l1);
                }
                else
                {
                    a = 0;
                    // std::cout << "Log Lik decreased";
                    // exit;
                }
                // myfile << "Aitken: " << a << std::endl;
                long double add_to = (1.0L / (1.0L - a)) * (l3 - l2);
                aloglik[it - 1] = loglik[it - 2] + add_to;
                //value below is arbitrary, had to be changed after modification in calculating Z_ig matrix and log likelihood
                if (abs_ld(aloglik[it - 1] - loglik[it - 2]) < 0.00000001L)
                {
                    exit_code = "Aitken's acceleration converged";
                    checks = 2;
                }
            }
        }
        // print it to R
        if (it % 10 == 0)
        {
            Rcpp::Rcout << it << std::endl;
        }
        if (it >= it_max)
        {
            exit_code = "Max iterations reached";
            checks = 999;
        }
        it++;
    }
    // myfile << "Exit code: " << checks << std::endl;
    // myfile.close();

    Rcpp::List result;
    result["loglik"] = Rcpp::wrap(loglik);
    result["z"] = Rcpp::wrap(z_loop);
    result["pi_g"] = Rcpp::wrap(pi_g);
    result["mu"] = Rcpp::wrap(mu);
    result["sigma"] = Rcpp::wrap(sigma);
    result["exit_code"] = Rcpp::wrap(exit_code);
    return Rcpp::wrap(result);
}
