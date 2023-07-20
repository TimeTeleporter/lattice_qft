#![feature(result_option_inspect)]
#![feature(drain_filter)]
#![feature(let_chains)]

use std::{error::Error, time::Instant};

use lattice_qft::{
    calculate_correlation_length,
    computation::ComputationSummary,
    export::{CorrelationLengths, CsvData, FitResult},
    kahan::WelfordsAlgorithm64,
};
use nalgebra::DVector;
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use rayon::prelude::*;
use varpro::{
    prelude::*,
    solvers::levmar::{LevMarProblemBuilder, LevMarSolver},
};

// Parameters
const BIN_SIZES: [usize; 6] = [1, 2, 4, 8, 16, 32];
const BOOTSTRAPPING_ENSEMBLES: u64 = 2000;
const SEEDED: bool = true;

// Analysis steps
const SYMMETRIZE: bool = true;
const CORR_FN_ANALYSIS: bool = false;
const FITTING: bool = false;
const COMP_CORR_FN_ANALYSIS: bool = true;

fn main() {
    let time = Instant::now();

    // Symmetrize the correlation functions
    if SYMMETRIZE {
        symmetrize_correlation_functions(
            lattice_qft::RESULTS_PATH,
            lattice_qft::CORR_FN_STATS_PATH_INCOMPLETE,
        );
    }

    // Binned resampling with compounding
    if CORR_FN_ANALYSIS {
        calculate_correlation_function_statistics_compounded(
            lattice_qft::RESULTS_PATH,
            lattice_qft::CORR_FN_STATS_PATH_INCOMPLETE,
            None,
        );
    }

    // Calculate the fits for each correlation function
    if FITTING {
        nonlin_fit(
            lattice_qft::RESULTS_PATH,
            lattice_qft::CORR_FN_STATS_PATH_INCOMPLETE,
        );
    }

    println!(
        "Correlation function analysis took {} secs",
        (time.elapsed().as_secs_f32() * 100.0).round() / 100.0
    );

    // Binned resampling with compounding
    if COMP_CORR_FN_ANALYSIS {
        calculate_correlation_function_statistics_compounded(
            lattice_qft::RESULTS_PATH,
            lattice_qft::COMP_CORR_FN_STATS_PATH_INCOMPLETE,
            Some(lattice_qft::COMP_RESULTS_PATH),
        );

        // Fit the compounded correlation functions
        nonlin_fit(
            lattice_qft::COMP_RESULTS_PATH,
            lattice_qft::COMP_CORR_FN_STATS_PATH_INCOMPLETE,
        );
    }

    println!(
        "Complete Analyis took {} secs",
        (time.elapsed().as_secs_f32() * 100.0).round() / 100.0
    );
}

/// This method calculates the statistics of the correlation functions of the simulations via binned bootstrapping resampling.
fn calculate_correlation_function_statistics_compounded(
    results_path: &str,
    comp_corr_fn_stats_path_incomplete: &str,
    comp_results_path: Option<&str>,
) {
    // Reading the results data
    let results = match ComputationSummary::fetch_csv_data(results_path, true) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Reading results for binned resampling: {}", err);
            return;
        }
    };

    let mut results: Vec<ComputationSummary> = results
        .into_iter()
        .filter(|summary| summary.correlation_data)
        .rev()
        .collect();

    let mut comp_results: Vec<ComputationSummary> = Vec::new();

    while let Some(mut key_summary) = results.pop() {
        key_summary.index;

        let mut local_results: Vec<ComputationSummary> = results
            .drain_filter(|summary| {
                comp_results_path.is_some()
                    && summary.temp == key_summary.temp
                    && summary.t == key_summary.t
                //&& summary.comptype == key_summary.comptype
            })
            .collect();

        local_results.reverse();
        local_results.push(key_summary.clone());
        local_results.reverse();

        let compund_data_amount: usize = local_results.len();

        // Fetching the correlation functions
        let local_corr_fns: Vec<Vec<Vec<f64>>> = local_results
            .into_iter()
            .filter_map(|summary| {
                if SYMMETRIZE {
                    fetch_correlation_functions_sym(summary.index)
                } else if CORR_FN_ANALYSIS {
                    fetch_correlation_functions(summary.index)
                } else {
                    Err("Unable to get correlation functions for analysis".into())
                }
                .inspect_err(|err| eprintln!("{}", err))
                .ok()
                .map(|corr_fns| corr_fns)
            })
            .collect();

        let bin_observations_iter = BIN_SIZES.into_iter().map(|bin_size| {
            // Binning the correlation functions
            let binned_local_corr_fns: Vec<(Vec<&[Vec<f64>]>, usize)> = local_corr_fns
                .iter()
                .map(|corr_fns| {
                    let binned_corr_fns: Vec<&[Vec<f64>]> =
                        corr_fns.chunks_exact(bin_size).collect();

                    // Quick check to see if we catch all correlation functions
                    let bin_samples_amount: usize = corr_fns.len() / bin_size;
                    assert_eq!(binned_corr_fns.len(), bin_samples_amount);

                    (binned_corr_fns, bin_samples_amount)
                })
                .collect();

            // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

            // We calculate the bootstrapping ensembles for each simulation and then combine them to extract the statistics.
            let bootstrapping_ensembles: Vec<(Vec<f64>, f64)> = (0..BOOTSTRAPPING_ENSEMBLES)
                .into_par_iter()
                .map(|seed| {
                    let mut rng: StdRng = if SEEDED {
                        StdRng::seed_from_u64(seed)
                    } else {
                        StdRng::from_entropy()
                    };

                    let binned_local_resampled_corr_fns: Vec<&(Vec<&[Vec<f64>]>, usize)> = (0
                        ..binned_local_corr_fns.len())
                        .into_iter()
                        .filter_map(|_| binned_local_corr_fns.choose(&mut rng))
                        .collect();

                    let single_ensemble_corr_fn_stats: Vec<WelfordsAlgorithm64> =
                        binned_local_resampled_corr_fns
                            .iter()
                            .map(|(binned_corr_fns, bin_samples_amount)| {
                                let single_sim_resampled_corr_fns: Vec<&Vec<f64>> = (0
                                    ..*bin_samples_amount)
                                    .into_iter()
                                    .filter_map(|_| binned_corr_fns.choose(&mut rng))
                                    .cloned()
                                    .flatten()
                                    .collect();
                                single_sim_resampled_corr_fns
                            })
                            .flatten()
                            .fold(
                                vec![WelfordsAlgorithm64::new(); key_summary.t],
                                |mut accumulator, single_sim_resampled_corr_fn| {
                                    accumulator
                                        .iter_mut()
                                        .zip(single_sim_resampled_corr_fn.into_iter())
                                        .for_each(|(accu, value)| accu.update(*value));
                                    accumulator
                                },
                            );

                    let single_ensemble_corr_fn_means: Vec<f64> = single_ensemble_corr_fn_stats
                        .into_iter()
                        .map(|corr_fn_value_stats| corr_fn_value_stats.get_mean())
                        .collect();

                    let p1: f64 = 2.0 * std::f64::consts::PI / (key_summary.t as f64);

                    let (single_ensemble_corr12, _): (f64, _) =
                        calculate_correlation_length(&single_ensemble_corr_fn_means, p1, p1 * 2.0);

                    (single_ensemble_corr_fn_means, single_ensemble_corr12)
                })
                .collect();

            let (corr_fn_stats, corr12_stats) = bootstrapping_ensembles.into_iter().fold(
                (
                    vec![WelfordsAlgorithm64::new(); key_summary.t],
                    WelfordsAlgorithm64::new(),
                ),
                |(mut corr_fn_accu, mut corr_len_accu), (corr_fn, corr_len)| {
                    corr_fn_accu
                        .iter_mut()
                        .zip(corr_fn.into_iter())
                        .for_each(|(accu, value)| accu.update(value));
                    corr_len_accu.update(corr_len);
                    (corr_fn_accu, corr_len_accu)
                },
            );

            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            let (corr_fn_means, corr_fn_sds) = corr_fn_stats
                .into_iter()
                .map(|stats| (stats.get_mean(), stats.get_sd_unbiased()))
                .unzip();

            let (corr12_mean, corr12_sd) =
                (corr12_stats.get_mean(), corr12_stats.get_sd_unbiased());

            (
                bin_size as u64,
                corr_fn_means,
                corr_fn_sds,
                corr12_mean,
                corr12_sd,
            )
        });

        let observations: Vec<(u64, Vec<f64>, Vec<f64>, f64, f64)> =
            bin_observations_iter.collect();

        // Filtering out the correlation function means
        let (corr_fn_means, remaining): (
            Vec<(u64, Vec<f64>)>,
            Vec<((u64, Vec<f64>), (u64, f64, f64))>,
        ) = observations
            .into_iter()
            .map(|resampled| {
                let (bin_size, corr_fn_means, corr_fn_sds, corr12_mean, corr12_sd) = resampled;
                let corr_fn_means = (bin_size, corr_fn_means);
                let corr_fn_sds = (bin_size, corr_fn_sds);
                let corr12 = (bin_size, corr12_mean, corr12_sd);
                (corr_fn_means, (corr_fn_sds, corr12))
            })
            .unzip();

        // Writing the correlation function means
        let path: &str = &(comp_corr_fn_stats_path_incomplete.to_owned()
            + &"correlation_"
            + &key_summary.index.to_string()
            + &"_res.csv");
        if let Err(err) = corr_fn_means.overwrite_csv(path) {
            eprintln!("Writing correlation function results: {}", err);
        }

        // Filtering out the correlation function errors
        let (corr_fn_sds, corr12s): (Vec<(u64, Vec<f64>)>, Vec<(u64, f64, f64)>) = remaining
            .into_iter()
            .map(|remaining| {
                let (corr_fn_sds, corr12) = remaining;
                (corr_fn_sds, corr12)
            })
            .unzip();

        // Writing the correlation fuction errors
        let path: &str = &(comp_corr_fn_stats_path_incomplete.to_owned()
            + &"correlation_"
            + &key_summary.index.to_string()
            + &"_err.csv");
        if let Err(err) = corr_fn_sds.overwrite_csv(path) {
            eprintln!("Writing correlation function errors: {}", err);
        }

        //Writing the correlation lenght data
        let correlation_lenghts: Vec<CorrelationLengths> = corr12s
            .into_iter()
            .map(|(bin_size, m12, m12_err)| {
                CorrelationLengths::new(key_summary.index, bin_size, m12).set_corr12_error(m12_err)
            })
            .collect();
        let path: &str = &(comp_corr_fn_stats_path_incomplete.to_owned()
            + &"correlation_"
            + &key_summary.index.to_string()
            + &"_len.csv");
        if let Err(err) = correlation_lenghts.overwrite_csv(path) {
            eprintln!("Writing correlation lengths: {}", err);
        }

        key_summary.comptype = compund_data_amount.to_string() + &" " + &key_summary.comptype + "s";

        comp_results.push(key_summary);
    }

    if let Some(path) = comp_results_path {
        if let Err(err) = comp_results.overwrite_csv(path) {
            eprintln!("Writing correlation lengths: {}", err);
        }
    }
}

/// Fetching the symmetrized correlation functions data for a given index
fn fetch_correlation_functions(index: u64) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
    let path: &str = &(lattice_qft::CORR_FN_PATH_INCOMPLETE.to_owned()
        + &"correlation_"
        + &index.to_string()
        + &".csv");
    let corr_fn: Result<Vec<Vec<f64>>, Box<dyn Error>> = Vec::<f64>::fetch_csv_data(path, false)
        .map_err(|err| format!("Fetching correlation functions {}: {}", path, err).into());
    corr_fn
}

fn symmetrize_correlation_functions(results_path: &str, corr_fn_stats_path_incomplete: &str) {
    // Reading the results data
    let results = match ComputationSummary::fetch_csv_data(results_path, true) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Reading results for symmetrizing: {}", err);
            return;
        }
    };

    let corr_fn_syms: Vec<(u64, Vec<Vec<f64>>)> = results
        .into_par_iter()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            // Reading the correlation functions
            let correlation_functions: Option<Vec<Vec<f64>>> =
                fetch_correlation_functions(summary.index)
                    .inspect_err(|err| eprintln!("{}", err))
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
        let path: &str = &(corr_fn_stats_path_incomplete.to_owned()
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

/// Fetching the symmetrized correlation functions data for a given index
fn fetch_correlation_functions_sym(index: u64) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
    let path: &str = &(lattice_qft::CORR_FN_STATS_PATH_INCOMPLETE.to_owned()
        + &"correlation_"
        + &index.to_string()
        + &"_sym.csv");
    let corr_fn: Result<Vec<Vec<f64>>, Box<dyn Error>> = Vec::<f64>::fetch_csv_data(path, false)
        .map_err(|err| format!("Fetching correlation functions {}: {}", path, err).into());
    corr_fn
}

/// Fetchin the correlation functions for a given index
fn fetch_correlation_functions_results(
    index: u64,
    corr_fn_stats_path_incomplete: &str,
) -> Result<Vec<(u64, Vec<f64>)>, Box<dyn Error>> {
    let path: &str = &(corr_fn_stats_path_incomplete.to_owned()
        + &"correlation_"
        + &index.to_string()
        + &"_res.csv");
    let corr_fns: Result<Vec<(u64, Vec<f64>)>, Box<dyn Error>> =
        <(u64, Vec<f64>) as CsvData>::fetch_csv_data(path, false)
            .map_err(|err| format!("Fetching {}: {}", path, err).into());
    corr_fns
}

fn nonlin_fit(results_path: &str, corr_fn_stats_path_incomplete: &str) {
    // Read the result data
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(results_path, true) {
            Ok(res) => res,
            Err(err) => {
                eprint!("Reading results for nonlinear fit: {}", err);
                return;
            }
        };

    let fitted: Vec<(u64, Vec<FitResult>)> = results
        .into_iter()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            fetch_correlation_functions_results(summary.index, corr_fn_stats_path_incomplete)
                .inspect_err(|err| eprintln!("Fitting: {}", err))
                .ok()
                .map(|x| (summary.index, x))
        })
        .map(|(index, corr_fn)| {
            let fits: Vec<FitResult> = corr_fn
                .into_iter()
                .filter_map(|(bin_size, corr_fn)| {
                    nonlin_regression(index, bin_size, corr_fn)
                        .map_err(|err| {
                            eprint!("{}", err);
                        })
                        //.or::<Box<dyn Error>>(Ok(FitResult::new(index, f64::NAN, f64::NAN, f64::NAN, f64::NAN)))
                        .ok()
                })
                .collect();
            (index, fits)
        })
        .collect();

    for (index, fits) in fitted {
        let path: &str = &(corr_fn_stats_path_incomplete.to_owned()
            + &"correlation_"
            + &index.to_string()
            + &"_fits.csv");
        if let Err(err) = fits.overwrite_csv(path) {
            eprintln!("Writing fits: {}", err);
        }
    }
}

fn nonlin_regression(
    index: u64,
    bin_size: u64,
    corr_fn: Vec<f64>,
) -> Result<FitResult, Box<dyn Error>> {
    match std::panic::catch_unwind(move || {
        let y_values: Vec<f64> = corr_fn;
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
        FitResult::new(index, bin_size, alpha[0], n, coeff[0], coeff[1])
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
