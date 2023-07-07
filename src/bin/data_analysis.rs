#![feature(result_option_inspect)]
#![feature(drain_filter)]
#![feature(let_chains)]

use std::error::Error;

use lattice_qft::{
    computation::ComputationSummary,
    export::{
        get_correlation_fn, get_correlation_fn_with_err, CorrelationLengths, CsvData, FitResult,
    },
    kahan::KahanSummation,
};
use nalgebra::DVector;
use varpro::{
    prelude::*,
    solvers::levmar::{LevMarProblemBuilder, LevMarSolver},
};

const SYMMETRIZE: bool = true;

fn main() {
    // Calculate the second moment correlation lenght for each correlation function
    second_moment(
        lattice_qft::RESULTS_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
        lattice_qft::RESULTS_SECOND_MOMENT_PATH,
        true,
    );
    // Calculate the fits for each correlation function
    nonlin_fit(
        lattice_qft::RESULTS_PATH,
        lattice_qft::CORR_FN_PATH_INCOMPLETE,
        lattice_qft::RESULTS_FIT_PATH,
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
    );
    // Fit the compounded correlation functions
    nonlin_fit(
        lattice_qft::RESULTS_COMP_PATH,
        lattice_qft::COMP_CORR_FN_PATH_INCOMPLETE,
        lattice_qft::RESULTS_COMP_FIT_PATH,
    );
}

fn compound_data(
    results_path: &str,
    results_comp_path: &str,
    corr_fn_path: &str,
    comp_corr_fn_path: &str,
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
    results.reverse();

    // Initialize new array for summaries
    let mut summaries: Vec<ComputationSummary> = Vec::new();

    // For each uniue set of coupling constant, lattice size and algorithm,
    // average the correlation functions and correlation lengths
    while let Some(summary) = results.pop() {
        let temp: Option<f64> = summary.temp;
        let max_t: Option<usize> = summary.t;
        let index: usize = summary.index;
        let comptype: Option<String> = summary.comptype.clone();

        // Read the correlation functions
        let mut corr_fn: Option<Vec<KahanSummation<f64>>> = get_correlation_fn(index, corr_fn_path)
            .map_err(|err| eprintln!("Index {} no correlation function: {}", index, err))
            .ok()
            .map(|ary| {
                ary.into_iter()
                    .map(|x| {
                        let mut kahan: KahanSummation<f64> = KahanSummation::new();
                        kahan.add(x);
                        kahan
                    })
                    .collect()
            });

        let mut corr: KahanSummation<f64> = KahanSummation::new();
        if let Some(corr12) = summary.corr12 {
            corr.add(corr12);
        } else {
            eprintln!("Index {} no correlation lenght availible!", index)
        }

        results
            .drain_filter(|entry| {
                entry.temp == temp && entry.t == max_t && entry.comptype == comptype
            })
            .for_each(|entry| {
                if let Some(corr12) = entry.corr12 {
                    corr.add(corr12);
                } else {
                    eprintln!("Index {} no correlation lenght availible!", entry.index)
                }
                if let Some(entry_corr_fn) = get_correlation_fn(entry.index, corr_fn_path)
                    .map_err(|err| {
                        eprintln!("Index {} no correlation function: {}", summary.index, err)
                    })
                    .ok() && let Some(corr_fn) = &mut corr_fn
                {
                    corr_fn
                        .iter_mut()
                        .zip(entry_corr_fn)
                        .for_each(|(corr1, corr2)| corr1.add(corr2))
                }
            });

        if let Some(corr_fn) = corr_fn {
            let corr_fn: Vec<f64> = corr_fn.into_iter().map(|kahan| kahan.mean()).collect();
            let path: &str =
                &(comp_corr_fn_path.to_owned() + &"correlation_" + &index.to_string() + &".csv");
            if let Err(err) = corr_fn.overwrite_csv(path) {
                eprintln!("Writing compounded correlation function: {}", err);
            }
        }

        summaries.push(summary.set_correlation_length(corr.mean()))
    }

    // Writing the compounded summary data
    if let Err(err) = summaries.clone().overwrite_csv(results_comp_path) {
        eprint!("Writing compounded summaries: {}", err);
    }

    // Building the errors for the compounded correlation functions
    if let Err(err) = build_corr_fn_errors(
        results_path,
        results_comp_path,
        corr_fn_path,
        comp_corr_fn_path,
    ) {
        eprintln!("Building correlation function errors: {}", err);
    }
}

fn build_corr_fn_errors(
    results_path: &str,
    results_comp_path: &str,
    corr_fn_path: &str,
    comp_corr_fn_path: &str,
) -> Result<(), Box<dyn Error>> {
    // Read the result data
    let mut results: Vec<ComputationSummary> =
        ComputationSummary::fetch_csv_data(results_path, true)?;
    let mut results_comp: Vec<ComputationSummary> =
        ComputationSummary::fetch_csv_data(results_comp_path, true)?;

    // For each compounded correlation function...
    while let Some(summary) = results_comp.pop() {
        let mean_comp_corr_fn: Vec<f64> = get_correlation_fn(summary.index, comp_corr_fn_path)?;
        let mut comp_corr_fn_errors: Vec<KahanSummation<f64>> = mean_comp_corr_fn
            .iter()
            .map(|_| KahanSummation::new())
            .collect();

        results
            .drain_filter(|entry| {
                entry.temp == summary.temp
                    && entry.t == summary.t
                    && entry.comptype == summary.comptype
            })
            .for_each(|entry| {
                let corr_fn: Vec<f64> = get_correlation_fn(entry.index, corr_fn_path).unwrap();
                comp_corr_fn_errors
                    .iter_mut()
                    .zip(corr_fn.into_iter().zip(mean_comp_corr_fn.iter()))
                    .for_each(|(kahan, (x, mean_x))| kahan.add((x - *mean_x) * (x - *mean_x)))
            });

        let comp_corr_fn_errors: Vec<f64> = comp_corr_fn_errors
            .into_iter()
            .map(|kahan| kahan.mean().sqrt())
            .collect();

        let path: &str = &(comp_corr_fn_path.to_owned()
            + &"correlation_"
            + &summary.index.to_string()
            + &"_err.csv");
        comp_corr_fn_errors.overwrite_csv(path)?;
    }

    Ok(())
}

fn nonlin_fit(results_path: &str, corr_fn_path: &str, fit_path: &str) {
    // Read the result data
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(results_path, true) {
            Ok(res) => res,
            Err(err) => {
                eprint!("{}", err);
                return;
            }
        };

    let fitted: Vec<FitResult> = results
        .into_iter()
        .filter(|res| res.correlation_data)
        .filter_map(|res| {
            nonlin_regression(res.index, corr_fn_path)
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

fn nonlin_regression(index: usize, corr_fn_path: &str) -> Result<FitResult, Box<dyn Error>> {
    match std::panic::catch_unwind(move || {
        let mut y_values = get_correlation_fn(index, corr_fn_path).unwrap();
        let n: usize = y_values.len();

        if SYMMETRIZE {
            let mut new_y_values: Vec<f64> = Vec::new();
            new_y_values.push(y_values[0]);
            for i in 1..n {
                new_y_values.push((y_values[i] + y_values[n - i]) / 2.0)
            }
            assert_eq!(n, new_y_values.len());
            y_values = new_y_values;
        }

        let n: f64 = n as f64;

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
            correlation_lenght_calculation(summary.index, summary.t?, corr_fn_path)
                .map_err(|err| {
                    eprint!("{}", err);
                })
                .ok()
                .inspect(|lenght| summary.corr12 = Some(lenght.corr12))
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
    index: usize,
    max_t: usize,
    incomplete_path: &str,
) -> Result<CorrelationLengths, Box<dyn Error>> {
    let (corr_fn, corr_fn_err): (Vec<f64>, Option<Vec<f64>>) =
        get_correlation_fn_with_err(index, incomplete_path)?;

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

    if let Some(corr_fn_err) = corr_fn_err {
        let (_, (m12_err, _)) =
            lattice_qft::calculate_correlation_length_errors(&corr_fn, &corr_fn_err, p1, 2.0 * p1);
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

    Ok(correlation_lengths)
}
