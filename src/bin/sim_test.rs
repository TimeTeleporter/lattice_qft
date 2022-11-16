#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use lattice_qft::{
    computation::{
        Algorithm, Computatable, Computation, ComputationResults, Observable, Simulation3d, Test3d,
    },
    export::CsvData,
    lattice::Lattice3d,
    REL_TEMP_ARY,
};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

const RANGE: usize = 16;

// Parameters for the Wilson loop
const _TEMP: f64 = 0.1;
const WIDTH: usize = 1;
const HEIGHT: usize = 1;
const BURNIN: usize = 100_000;
const ITERATIONS: usize = 10_000_000;

const RESULTS_PATH: &str = "./data/results.csv";

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    // Initialise the simulations
    let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
    for temp in REL_TEMP_ARY {
        let observable: Observable = Observable::Wilson {
            width: WIDTH,
            height: HEIGHT,
        };
        //comps.push(Test3d::new_computation(&lattice, observable, temp, RANGE));
        comps.push(Simulation3d::new_compuatation(
            &lattice,
            Algorithm::Cluster,
            observable,
            temp,
            BURNIN,
            ITERATIONS,
        ));
        comps.push(Simulation3d::new_compuatation(
            &lattice,
            Algorithm::Metropolis,
            observable,
            temp,
            BURNIN,
            ITERATIONS,
        ));
    }

    // Parallel over all temp data
    let data: Vec<ComputationResults> = comps
        .into_par_iter()
        .filter_map(|comp| comp.run().ok())
        .collect();

    println!("All calcualtions done");

    for entry in data {
        if let Err(err) = entry.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        };
    }
}
