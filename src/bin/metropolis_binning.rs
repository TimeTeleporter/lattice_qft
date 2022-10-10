//! Performs a Metropolis simulation for a given temperature. The resulting
//! Markov chain gets binned repeatedly to eliminate the correlation between
//! the calculated observables in order calculate the real error of the
//! observable.
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use std::error::Error;

use serde::{Deserialize, Serialize};

use lattice_qft::{
    export::{clean_csv, CsvData},
    lattice::Lattice3d,
    metropolis::{metropolis_simulation3d, MarkovChain, MetropolisSimResult},
};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;

const BINS_PATH: &str = "data/metropolis_binning/bins1_58.csv";
const RESULTS_PATH: &str = "data/metropolis_binning/binning_sim_results.csv";

const BURNIN: usize = 100_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 10_000_000;

const TEMP: f64 = 1.58; // Make sure that the path corresponds to the temp

/// Datatype to save and read binned simulation output.
#[derive(Debug, Serialize, Deserialize)]
struct MetropolisBinData {
    temp: f64,
    binning_step: u32,
    variance: f64, // The square of the standard deviation
    std: f64,
    error: f64,
}

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    let (result, data) = metropolis_simulation3d(&lattice, TEMP, BURNIN, ITERATIONS);

    if let Err(err) = metropolis_output(data, result, TEMP) {
        eprint!("{}", err);
    }

    println!("{}: Done", TEMP);
}

fn metropolis_output(
    data: MarkovChain,
    result: MetropolisSimResult,
    temp: f64,
) -> Result<(), Box<dyn Error>> {
    result.clone().read_write_csv(RESULTS_PATH)?;

    clean_csv(BINS_PATH)?;

    for bin_variance in calculate_bin_var(data, temp) {
        bin_variance.read_write_csv(BINS_PATH)?;
    }

    Ok(())
}

fn calculate_bin_var(mut data: MarkovChain, temp: f64) -> Vec<MetropolisBinData> {
    let mut results: Vec<MetropolisBinData> = Vec::new();
    let original_length: usize = data.len();

    let mut binning_step: u32 = 0;
    while 2_usize.pow(binning_step) < original_length {
        // Calculating the relevant values
        let variance: f64 = data.calulate_variance();
        let std: f64 = variance.sqrt();
        let error: f64 = std / ((data.len() - 1) as f64).sqrt();

        // Printing
        results.push(MetropolisBinData {
            temp,
            binning_step,
            variance,
            std,
            error,
        });

        // Binning
        data = data.binning();
        binning_step += 1;
    }

    results
}
