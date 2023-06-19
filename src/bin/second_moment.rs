use std::error::Error;

use lattice_qft::{
    computation::ComputationSummary,
    export::{clean_csv, get_correlation_fn, CorrelationLengths, CsvData},
};

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

    if let Err(err) = clean_csv(lattice_qft::RESULTS_CORR_PATH)
        .and(fitted.read_write_csv(lattice_qft::RESULTS_CORR_PATH, true))
    {
        eprint!("{}", err);
    }
}

fn correlation_lenght_calculation(
    index: usize,
    max_t: usize,
) -> Result<CorrelationLengths, Box<dyn Error>> {
    let corr_fn: Vec<f64> = get_correlation_fn(index)?;

    let (m12, m23, m13) = return_correlation_lengths(corr_fn, max_t)?;

    Ok(CorrelationLengths::new(index, m12, m23, m13))
}

fn return_correlation_lengths(
    correlation_fn: Vec<f64>,
    max_t: usize,
) -> Result<(f64, f64, f64), Box<dyn Error>> {
    let Some(x_max): Option<f64> = correlation_fn.clone().into_iter().reduce(f64::max) else {
        return Err("Unable to determine the  maximum of the correlation function".into());
    };

    let correlation_fn: Vec<f64> = correlation_fn.into_iter().map(|x| x_max - x).collect();

    let p1: f64 = 2.0 * std::f64::consts::PI / (max_t as f64);

    let m12: f64 = calculate_correlation_length(correlation_fn.clone(), p1, 2.0 * p1);
    let m23: f64 = calculate_correlation_length(correlation_fn.clone(), 2.0 * p1, 3.0 * p1);
    let m13: f64 = calculate_correlation_length(correlation_fn, p1, 3.0 * p1);

    Ok((m12, m23, m13))
}

fn calculate_correlation_length(correlation_fn: Vec<f64>, p1: f64, p2: f64) -> f64 {
    let (g1_re, g1_im): (f64, f64) = discrete_fourier_transform(correlation_fn.clone(), p1);
    let (g2_re, g2_im): (f64, f64) = discrete_fourier_transform(correlation_fn, p2);

    let cos1: f64 = p1.cos();
    let cos2: f64 = p2.cos();

    let g1_re2: f64 = g1_re * g1_re;
    let g2_re2: f64 = g2_re * g2_re;
    let g1_im2: f64 = g1_im * g1_im;
    let g2_im2: f64 = g2_im * g2_im;
    let g_re_mix: f64 = g1_re * g2_re;
    let g_im_mix: f64 = g1_im * g2_im;

    let abs: f64 = g1_re2 - 2.0 * g_re_mix + g2_re2 + g1_im2 - 2.0 * g_im_mix + g2_im2;

    let cosh_re: f64 =
        g1_re2 * cos1 - g_re_mix * cos1 - g_re_mix * cos2 + g2_re2 * cos2 + g1_im2 * cos1
            - g_im_mix * cos1
            - g_im_mix * cos2
            + g2_im2 * cos2;

    let cosh = cosh_re / abs;

    assert!(!cosh.is_nan());

    f64::acosh(cosh)
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
