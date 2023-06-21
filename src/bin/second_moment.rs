use std::error::Error;

use lattice_qft::{
    computation::ComputationSummary,
    export::{get_correlation_fn, CorrelationLengths, CsvData},
};
use num_complex::Complex64;

fn main() {
    // Read the result data
    let results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true) {
            Ok(summary) => summary,
            Err(err) => {
                eprint!("{}", err);
                return;
            }
        };

    let fitted: Vec<CorrelationLengths> = results
        .into_iter()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            correlation_lenght_calculation(summary.index, summary.t?)
                .map_err(|err| {
                    eprint!("{}", err);
                })
                .ok()
        })
        .collect();

    if let Err(err) = fitted.overwrite_csv(lattice_qft::RESULTS_CORR_PATH) {
        eprint!("{}", err);
    }
}

fn correlation_lenght_calculation(
    index: usize,
    max_t: usize,
) -> Result<CorrelationLengths, Box<dyn Error>> {
    let corr_fn: Vec<f64> = get_correlation_fn(index)?;

    let ary = return_correlation_lengths(corr_fn, max_t)?;

    Ok(CorrelationLengths::new(index, ary))
}

fn return_correlation_lengths(
    correlation_fn: Vec<f64>,
    max_t: usize,
) -> Result<[f64; 6], Box<dyn Error>> {
    let Some(x_max): Option<f64> = correlation_fn.clone().into_iter().reduce(f64::max) else {
        return Err("Unable to determine the  maximum of the correlation function".into());
    };

    let correlation_fn: Vec<f64> = correlation_fn.into_iter().map(|x| x_max - x).collect();

    let p1: f64 = 2.0 * std::f64::consts::PI / (max_t as f64);

    let m12: f64 = calculate_correlation_length(correlation_fn.clone(), p1, 2.0 * p1).0;
    let m23: f64 = calculate_correlation_length(correlation_fn.clone(), 2.0 * p1, 3.0 * p1).0;
    let m34: f64 = calculate_correlation_length(correlation_fn.clone(), 3.0 * p1, 4.0 * p1).0;
    let m13: f64 = calculate_correlation_length(correlation_fn.clone(), 1.0 * p1, 3.0 * p1).0;
    let m24: f64 = calculate_correlation_length(correlation_fn.clone(), 2.0 * p1, 4.0 * p1).0;
    let m14: f64 = calculate_correlation_length(correlation_fn.clone(), 1.0 * p1, 4.0 * p1).0;

    Ok([m12, m23, m34, m13, m24, m14])
}

fn calculate_correlation_length(correlation_fn: Vec<f64>, p1: f64, p2: f64) -> (f64, f64) {
    let (g1_re, g1_im): (f64, f64) = discrete_fourier_transform(correlation_fn.clone(), p1);
    let (g2_re, g2_im): (f64, f64) = discrete_fourier_transform(correlation_fn, p2);

    let g1: Complex64 = Complex64::new(g1_re, g1_im);
    let g2: Complex64 = Complex64::new(g2_re, g2_im);

    let cos1: f64 = p1.cos();
    let cos2: f64 = p2.cos();

    let cosh: Complex64 = (g1 * cos1 - g2 * cos2).fdiv(g1 - g2);

    assert!(!cosh.is_nan());

    let ma: Complex64 = cosh.acosh();

    (ma.re, ma.im)
}

fn discrete_fourier_transform(correlation_fn: Vec<f64>, momentum: f64) -> (f64, f64) {
    let real: f64 = correlation_fn
        .clone()
        .into_iter()
        .enumerate()
        .map(|(x, g_x)| g_x * f64::cos(momentum * (x as f64)))
        .sum();
    let imaginary: f64 = correlation_fn
        .into_iter()
        .enumerate()
        .map(|(x, g_x)| g_x * f64::sin(momentum * (x as f64)))
        .sum();
    (real, imaginary)
}
