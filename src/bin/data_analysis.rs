#![feature(result_option_inspect)]
#![feature(drain_filter)]
#![feature(let_chains)]

use std::error::Error;

use lattice_qft::{
    computation::ComputationSummary,
    export::{CorrelationLengths, CsvData, FitResult},
    kahan::WelfordsAlgorithm64,
};
use nalgebra::DVector;
use rand::{rngs::ThreadRng, seq::SliceRandom};
use rayon::prelude::*;
use varpro::{
    prelude::*,
    solvers::levmar::{LevMarProblemBuilder, LevMarSolver},
};

// Global constants
const MAX_T: usize = 54;

// Parameters
const BOOTSTRAPPING_SAMPLES: u64 = 100;

fn main() {
    // Symmetrize the correlation functions
    symmetrize_correlation_functions(
        lattice_qft::RESULTS_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
    );

    // Binned resampling of the correlation functions
    calculate_correlation_function_statistics(
        lattice_qft::RESULTS_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
        BOOTSTRAPPING_SAMPLES,
    );

    // Calculate the fits for each correlation function
    nonlin_fit(
        lattice_qft::RESULTS_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
        lattice_qft::CORR_FN_FIT_PATH,
    );

    /* Old data analysis code
    // Calculate the second moment correlation lenght for each correlation function
    second_moment(
        lattice_qft::RESULTS_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
        lattice_qft::RESULTS_SECOND_MOMENT_PATH,
        true,
        false,
    );

    // Calculated the compounded second moment correlation lenghts and correlation functions
    compound_data(
        lattice_qft::RESULTS_PATH,
        lattice_qft::RESULTS_COMP_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
        lattice_qft::COMP_CORR_FN_PATH_INCOMPLETE,
    );
    // Calculate the second moment correlation lenght from the compounded correlation functions
    second_moment(
        lattice_qft::RESULTS_COMP_PATH,
        lattice_qft::COMP_CORR_FN_PATH_INCOMPLETE,
        lattice_qft::RESULTS_COMP_SECOND_MOMENT_PATH,
        false,
        true,
    );
    // Fit the compounded correlation functions
    nonlin_fit(
        lattice_qft::RESULTS_COMP_PATH,
        lattice_qft::COMP_CORR_FN_PATH_INCOMPLETE,
        lattice_qft::RESULTS_COMP_FIT_PATH,
    );*/
}

fn symmetrize_correlation_functions(results_path: &str, corr_fn_path_incomplete: &str) {
    // Reading the results data
    let results = match ComputationSummary::fetch_csv_data(results_path, true) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Reading results for symmetrizing: {}", err);
            return;
        }
    };

    let corr_fn_syms: Vec<(u64, Vec<Vec<f64>>)> = results
        .into_iter()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            // Reading the correlation functions
            let path: &str = &(corr_fn_path_incomplete.to_owned()
                + &"correlation_"
                + &summary.index.to_string()
                + &".csv");
            let correlation_functions: Option<Vec<Vec<f64>>> =
                Vec::<f64>::fetch_csv_data(path, false)
                    .inspect_err(|err| {
                        eprintln!("Fetching correlation functions {}: {}", path, err)
                    })
                    .ok();
            correlation_functions
                .map(|correlation_functions| (summary.index, correlation_functions))
        })
        .map(|(index, mut correlation_functions)| {
            // Symmetrizing the correlation functions
            correlation_functions.iter_mut().for_each(|corr_fn| {
                symmetrize_single_correlation_function(corr_fn);
            });
            (index, correlation_functions)
        })
        .collect();

    // Writing the symmetrized correlation functions
    for (index, corr_fn_sym) in corr_fn_syms {
        let path: &str = &(corr_fn_path_incomplete.to_owned()
            + &"correlation_"
            + &index.to_string()
            + &"_sym.csv");
        if let Err(err) = corr_fn_sym.overwrite_csv(path) {
            eprintln!("Writing the symmetrized correlation function: {}", err);
        }
    }
}

fn symmetrize_single_correlation_function(corr_fn: &mut Vec<f64>) {
    let copied = corr_fn.clone();
    corr_fn
        .iter_mut()
        .skip(1)
        .zip(copied.into_iter().rev())
        .for_each(|(up, down)| *up = (*up + down) / 2.0);
}

#[test]
fn test_symmetrize_single_correlation_function() {
    let mut corr_fn: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let corr_fn_sym: Vec<f64> = vec![1.0, 4.0, 4.0, 4.0, 4.0, 4.0];
    symmetrize_single_correlation_function(&mut corr_fn);
    assert_eq!(corr_fn, corr_fn_sym);
}

/// Fetchin the symmetrized correlation functions data for a given index
fn fetch_correlation_functions_sym(
    index: u64,
    corr_fn_path_incomplete: &str,
) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
    let path: &str =
        &(corr_fn_path_incomplete.to_owned() + &"correlation_" + &index.to_string() + &"_sym.csv");
    let corr_fn: Result<Vec<Vec<f64>>, Box<dyn Error>> = Vec::<f64>::fetch_csv_data(path, false)
        .map_err(|err| format!("Fetching correlation functions {}: {}", path, err).into());
    corr_fn
}

/// This method calculates the statistics of the correlation functions of the simulations via binned bootstrapping resampling.
fn calculate_correlation_function_statistics(
    results_path: &str,
    corr_fn_path_incomplete: &str,
    bootstrapping_ensembles: u64,
) {
    // Reading the results data
    let results = match ComputationSummary::fetch_csv_data(results_path, true) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Reading results for binned resampling: {}", err);
            return;
        }
    };

    let bin_size: usize = 4;

    let (corr_fn_results, corr_fn_errors): (Vec<(u64, Vec<f64>)>, Vec<(u64, Vec<f64>)>) = results
        .par_iter()
        .filter(|summary| summary.correlation_data)
        // Fetching the correlation functions
        .filter_map(|summary| {
            fetch_correlation_functions_sym(summary.index, corr_fn_path_incomplete)
                .inspect_err(|err| eprintln!("{}", err))
                .ok()
                .map(|corr_fns| (summary, corr_fns))
        })
        // Next we create a bootstrap dataset
        .map(|(summary, corr_fns)| {
            let mut rng = ThreadRng::default();

            // Checking if the correlation function lenght corresponds to the lenght in the summary.
            assert_eq!(summary.t, corr_fns[0].len());

            let binned_corr_fns: Vec<&[Vec<f64>]> = corr_fns.chunks_exact(bin_size).collect();

            // Quick check to see if we catch all correlation functions
            let bin_samples_amount: usize = corr_fns.len() / bin_size;
            assert_eq!(binned_corr_fns.len(), bin_samples_amount);

            // Setting up the statistics collection
            let mut statistics: Vec<WelfordsAlgorithm64> = Vec::with_capacity(MAX_T);
            for _ in 0..summary.t {
                statistics.push(WelfordsAlgorithm64::new());
            }

            // We repeat the resampling to build new ensembles by...
            (0..bootstrapping_ensembles).into_iter().for_each(|_| {
                let mut sample_statistics: Vec<WelfordsAlgorithm64> = Vec::with_capacity(MAX_T);
                for _ in 0..summary.t {
                    sample_statistics.push(WelfordsAlgorithm64::new());
                }

                // ...drawing a random sample with replacement for each original sample. We employ
                // bins to combat the autocorrelation of the data.
                (0..bin_samples_amount)
                    .into_iter()
                    .filter_map(|_| binned_corr_fns.choose(&mut rng))
                    .cloned()
                    .flatten()
                    .for_each(|corr_fn| {
                        sample_statistics
                            .iter_mut()
                            .zip(corr_fn.iter())
                            .for_each(|(sample, value)| sample.update(*value))
                    });

                sample_statistics
                    .into_iter()
                    .map(|sample| sample.get_mean())
                    .zip(statistics.iter_mut())
                    .for_each(|(value, stat)| stat.update(value));
            });

            // For each ensemble we then do our analysis to get the real variance of the data.
            let (means, sds): (Vec<f64>, Vec<f64>) = statistics
                .iter()
                .map(|stat| (stat.get_mean(), stat.get_sd_unbiased()))
                .unzip();

            ((summary.index, means), (summary.index, sds))
        })
        .collect::<Vec<((u64, Vec<f64>), (u64, Vec<f64>))>>()
        .into_iter()
        .unzip();

    // Writing the mean correlation function
    for (index, corr_fn_result) in corr_fn_results {
        let path: &str = &(corr_fn_path_incomplete.to_owned()
            + &"correlation_"
            + &index.to_string()
            + &"_res.csv");
        if let Err(err) = corr_fn_result.overwrite_csv(path) {
            eprintln!("Writing correlation function results: {}", err);
        }
    }

    // Writing the correlation function err
    for (index, corr_fn_error) in corr_fn_errors {
        let path: &str = &(corr_fn_path_incomplete.to_owned()
            + &"correlation_"
            + &index.to_string()
            + &"_err.csv");
        if let Err(err) = corr_fn_error.overwrite_csv(path) {
            eprintln!("Writing correlation function errors: {}", err);
        }
    }
}

/// Fetchin the correlation functions for a given index
fn fetch_correlation_functions_results(
    index: u64,
    corr_fn_path_incomplete: &str,
) -> Result<Vec<f64>, Box<dyn Error>> {
    let path: &str =
        &(corr_fn_path_incomplete.to_owned() + &"correlation_" + &index.to_string() + &"_res.csv");
    let corr_fn: Result<Vec<f64>, Box<dyn Error>> = f64::fetch_csv_data(path, false)
        .map_err(|err| format!("Fetching {}: {}", path, err).into());
    corr_fn
}

/// Fetchin the errors of the correlation functions for a given index
fn fetch_correlation_functions_errors(
    index: u64,
    corr_fn_path_incomplete: &str,
) -> Result<Vec<f64>, Box<dyn Error>> {
    let path: &str =
        &(corr_fn_path_incomplete.to_owned() + &"correlation_" + &index.to_string() + &"_err.csv");
    let corr_fn: Result<Vec<f64>, Box<dyn Error>> = f64::fetch_csv_data(path, false)
        .map_err(|err| format!("Fetching {}: {}", path, err).into());
    corr_fn
}

fn nonlin_fit(results_path: &str, corr_fn_path_incomplete: &str, fit_path: &str) {
    // Read the result data
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(results_path, true) {
            Ok(res) => res,
            Err(err) => {
                eprint!("Reading results for nonlinear fit: {}", err);
                return;
            }
        };

    let fitted: Vec<FitResult> = results
        .into_iter()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            fetch_correlation_functions_results(summary.index, corr_fn_path_incomplete)
                .inspect_err(|err| eprintln!("Fitting: {}", err))
                .ok()
                .map(|x| (summary.index, x))
        })
        .filter_map(|(index, corr_fn)| {
            nonlin_regression(index, corr_fn)
                .map_err(|err| {
                    eprint!("{}", err);
                })
                //.or::<Box<dyn Error>>(Ok(FitResult::new(index, f64::NAN, f64::NAN, f64::NAN, f64::NAN)))
                .ok()
        })
        .collect();

    if let Err(err) = fitted.overwrite_csv(fit_path) {
        eprint!("{}", err);
    }
}

fn nonlin_regression(index: u64, corr_fn: Vec<f64>) -> Result<FitResult, Box<dyn Error>> {
    match std::panic::catch_unwind(move || {
        let mut y_values: Vec<f64> = corr_fn;
        let n: f64 = y_values.len() as f64;

        // Edit the data
        let Some(x_max): Option<f64> = y_values.clone().into_iter().reduce(f64::max) else {
            panic!("Unable to find the maximum");
        };

        let y_values: Vec<f64> = y_values.into_iter().map(|x| x_max - x).collect();

        let x_values: Vec<f64> = y_values.iter().enumerate().map(|(x, _)| x as f64).collect();

        let x = DVector::from_vec(x_values);
        let y = DVector::from_vec(y_values);

        // The function to fit
        fn nonlin_fn(x: &DVector<f64>, m: f64, n: f64) -> DVector<f64> {
            x.map(|x| (m * (x - n / 2.0)).cosh())
        }

        let nonlin_fn_n_given = move |x: &DVector<f64>, m: f64| nonlin_fn(x, m, n);

        // Derivative with respect to m
        fn nonlin_fn_dm(x: &DVector<f64>, m: f64, n: f64) -> DVector<f64> {
            x.map(|x| (m * (x - n / 2.0)).sinh() * (x - n / 2.0))
        }

        let nonlin_fn_dm_n_given = move |x: &DVector<f64>, m: f64| nonlin_fn_dm(x, m, n);

        /* Derivative with respect to n
        fn nonlin_fn_dn(x: &DVector<f64>, m: f64, n: f64) -> DVector<f64> {
            x.map(|x| -1.0 * (m * (x - n / 2.0)).sinh() * (m / 2.0))
        }*/

        let model = SeparableModelBuilder::<f64>::new(&["m"])
            .function(&["m"], nonlin_fn_n_given)
            .partial_deriv("m", nonlin_fn_dm_n_given)
            .invariant_function(|x| DVector::from_element(x.len(), 1.0))
            .independent_variable(x)
            .initial_parameters(vec![0.1])
            .build()
            .unwrap();

        let problem = LevMarProblemBuilder::new(model)
            .observations(y)
            .build()
            .unwrap();

        // This thing panics, we need to catch it.
        let (solved_problem, report) = LevMarSolver::new()
            .with_xtol(f64::EPSILON)
            .minimize(problem);
        if !report.termination.was_successful() {
            panic!("termination was not successful");
        }

        let alpha = solved_problem.params();
        let coeff = solved_problem
            .linear_coefficients()
            .ok_or("Unable to compute the coefficients")
            .unwrap();

        // FitResult::new(index, alpha[0], alpha[1], coeff[0], coeff[1])
        FitResult::new(index, alpha[0], n, coeff[0], coeff[1])
    }) {
        Ok(fit) => Ok(fit),
        Err(err) => {
            let err_msg = match (err.downcast_ref(), err.downcast_ref::<String>()) {
                (Some(&s), _) => s,
                (_, Some(s)) => &**s,
                _ => "<No panic message>",
            };
            Err(format!("nonlin regression paniced: {}", err_msg).into())
        }
    }
}

fn second_moment(
    results_path: &str,
    corr_fn_path: &str,
    second_moment_path: &str,
    write_corr12_in_results: bool,
    do_error_calculations: bool,
) {
    // Read the result data
    let mut results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(results_path, true) {
            Ok(summary) => summary,
            Err(err) => {
                eprint!("{}", err);
                return;
            }
        };

    // Calculate the second moment parameters
    let fitted: Vec<CorrelationLengths> = results
        .iter_mut()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            correlation_lenght_calculation(
                summary.index,
                summary.t,
                corr_fn_path,
                do_error_calculations,
            )
            .map_err(|err| {
                eprint!("{}", err);
            })
            .ok()
        })
        .collect();

    if let Err(err) = fitted.overwrite_csv(second_moment_path) {
        eprint!("{}", err);
    }
    if write_corr12_in_results {
        if let Err(err) = results.overwrite_csv(results_path) {
            eprint!("{}", err);
        }
    }
}

fn correlation_lenght_calculation(
    index: u64,
    max_t: usize,
    incomplete_path: &str,
    do_error_calculations: bool,
) -> Result<CorrelationLengths, Box<dyn Error>> {
    let (corr_fn, corr_fn_err): (Vec<f64>, Option<Vec<f64>>) = { todo!() };

    let Some(x_max): Option<f64> = corr_fn.clone().into_iter().reduce(f64::max) else {
            return Err("Unable to determine the  maximum of the correlation function".into());
        };

    let corr_fn: Vec<f64> = corr_fn.into_iter().map(|x| x_max - x).collect();
    let p1: f64 = 2.0 * std::f64::consts::PI / (max_t as f64);

    let m12: f64 = lattice_qft::calculate_correlation_length(&corr_fn, p1, 2.0 * p1).0;
    let m23: f64 = lattice_qft::calculate_correlation_length(&corr_fn, 2.0 * p1, 3.0 * p1).0;
    let m34: f64 = lattice_qft::calculate_correlation_length(&corr_fn, 3.0 * p1, 4.0 * p1).0;
    let m13: f64 = lattice_qft::calculate_correlation_length(&corr_fn, 1.0 * p1, 3.0 * p1).0;
    let m24: f64 = lattice_qft::calculate_correlation_length(&corr_fn, 2.0 * p1, 4.0 * p1).0;
    let m14: f64 = lattice_qft::calculate_correlation_length(&corr_fn, 1.0 * p1, 4.0 * p1).0;

    let values: [f64; 6] = [m12, m23, m34, m13, m24, m14];
    let mut correlation_lengths: CorrelationLengths = CorrelationLengths::new(index, values);

    if do_error_calculations {
        if let Some(corr_fn_err) = corr_fn_err {
            let (_, (m12_err, _)) = lattice_qft::calculate_correlation_length_errors(
                &corr_fn,
                &corr_fn_err,
                p1,
                2.0 * p1,
            );
            let (_, (m23_err, _)) = lattice_qft::calculate_correlation_length_errors(
                &corr_fn,
                &corr_fn_err,
                2.0 * p1,
                3.0 * p1,
            );
            let (_, (m34_err, _)) = lattice_qft::calculate_correlation_length_errors(
                &corr_fn,
                &corr_fn_err,
                3.0 * p1,
                4.0 * p1,
            );
            let (_, (m13_err, _)) = lattice_qft::calculate_correlation_length_errors(
                &corr_fn,
                &corr_fn_err,
                1.0 * p1,
                3.0 * p1,
            );
            let (_, (m24_err, _)) = lattice_qft::calculate_correlation_length_errors(
                &corr_fn,
                &corr_fn_err,
                2.0 * p1,
                4.0 * p1,
            );
            let (_, (m14_err, _)) = lattice_qft::calculate_correlation_length_errors(
                &corr_fn,
                &corr_fn_err,
                1.0 * p1,
                4.0 * p1,
            );

            let errors: [Option<f64>; 6] = [
                Some(m12_err),
                Some(m23_err),
                Some(m34_err),
                Some(m13_err),
                Some(m24_err),
                Some(m14_err),
            ];
            correlation_lengths.set_errors(errors);
        }
    }

    Ok(correlation_lengths)
}
