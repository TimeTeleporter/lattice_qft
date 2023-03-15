#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use lattice_qft::{
    //REL_TEMP_ARY,
    algorithm::AlgorithmType,
    computation::{Computation, ComputationResult, Compute},
    export::CsvData,
    lattice::Lattice3d,
    observable::ObservableType,
};

const TEST_X: usize = 4;
const TEST_Y: usize = 5;
const TEST_T: usize = 6;
const SIZE: usize = TEST_X * TEST_Y * TEST_T;

const RANGE: usize = 12;

const BURNIN: usize = 10_000;
const ITERATIONS: usize = 1_000_000;

const WIDTH: usize = TEST_X / 3;
const HEIGHT: usize = TEST_T;

const RESULTS_PATH: &str = "./data/results.csv";

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    // Initialise the simulations
    let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
    for temp in [0.1] {
        let observable: ObservableType = ObservableType::SizeNormalizedAction;
        comps.push(Computation::new_wilson_test(
            &lattice,
            temp,
            RANGE,
            WIDTH,
            HEIGHT,
            observable.clone(),
        ));
        comps.push(Computation::new_wilson_sim(
            &lattice,
            temp,
            AlgorithmType::new_metropolis(),
            BURNIN,
            ITERATIONS,
            WIDTH,
            HEIGHT,
            observable.clone(),
        ));
        comps.push(Computation::new_test(
            &lattice,
            temp,
            RANGE,
            observable.clone(),
        ));
        comps.push(Computation::new_simulation(
            &lattice,
            temp,
            AlgorithmType::new_metropolis(),
            BURNIN,
            ITERATIONS,
            observable,
        ));
    }

    // Parallel over all temp data
    let data: Vec<ComputationResult<3, SIZE>> = comps
        .into_par_iter()
        .filter_map(|comp| comp.run().ok())
        .collect();

    println!("All calcualtions done");

    for entry in data {
        if let Err(err) = entry.into_export().read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        };
    }
}

#[test]
fn test_const() {
    assert_eq!(SIZE, TEST_X * TEST_Y * TEST_T);
}
