#![feature(result_option_inspect)]

use std::error::Error;

use lattice_qft::{
    computation::ComputationSummary,
    export::{get_correlation_fn, CorrelationLengths, CsvData},
};

fn main() {
    // Read the result data
    let mut results: Vec<ComputationSummary> =
        match ComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true) {
            Ok(summary) => summary,
            Err(err) => {
                eprint!("{}", err);
                return;
            }
        };

    let fitted: Vec<CorrelationLengths> = results
        .iter_mut()
        .filter(|summary| summary.correlation_data)
        .filter_map(|summary| {
            correlation_lenght_calculation(
                summary.index,
                summary.t?,
                lattice_qft::CORRELATION_PLOT_PATH_INCOMPLETE,
            )
            .map_err(|err| {
                eprint!("{}", err);
            })
            .ok()
            .inspect(|lenght| summary.corr12 = Some(lenght.corr12))
        })
        .collect();

    if let Err(err) = fitted.overwrite_csv(lattice_qft::RESULTS_CORR_PATH) {
        eprint!("{}", err);
    }
    if let Err(err) = results.overwrite_csv(lattice_qft::RESULTS_PATH) {
        eprint!("{}", err);
    }
}

fn correlation_lenght_calculation(
    index: usize,
    max_t: usize,
    incomplete_path: &str,
) -> Result<CorrelationLengths, Box<dyn Error>> {
    let corr_fn: Vec<f64> = get_correlation_fn(index, incomplete_path)?;

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

    let m12: f64 =
        lattice_qft::outputdata::calculate_correlation_length(&correlation_fn, p1, 2.0 * p1).0;
    let m23: f64 =
        lattice_qft::outputdata::calculate_correlation_length(&correlation_fn, 2.0 * p1, 3.0 * p1)
            .0;
    let m34: f64 =
        lattice_qft::outputdata::calculate_correlation_length(&correlation_fn, 3.0 * p1, 4.0 * p1)
            .0;
    let m13: f64 =
        lattice_qft::outputdata::calculate_correlation_length(&correlation_fn, 1.0 * p1, 3.0 * p1)
            .0;
    let m24: f64 =
        lattice_qft::outputdata::calculate_correlation_length(&correlation_fn, 2.0 * p1, 4.0 * p1)
            .0;
    let m14: f64 =
        lattice_qft::outputdata::calculate_correlation_length(&correlation_fn, 1.0 * p1, 4.0 * p1)
            .0;

    Ok([m12, m23, m34, m13, m24, m14])
}
