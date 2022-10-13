//! This binary calculates a metropolis simulation and saves the result in a csv file.
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use lattice_qft::{
    export::CsvData, export::SimResult, lattice::Lattice3d, metropolis::metropolis_simulation3d,
    LONG_TEMP_ARY,
};

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;

const RESULTS_PATH: &str = "data/metropolis_test/sim_data_long.csv";

const BURNIN: usize = 100_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 10_000_000;

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    let results: Vec<SimResult> = LONG_TEMP_ARY
        .par_iter()
        .map(|&temp| {
            let (result, _) = metropolis_simulation3d(&lattice, temp, BURNIN, ITERATIONS);
            result
        })
        .collect();

    for res in results {
        if let Err(err) = res.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        }
    }

    println!("Done");
}
