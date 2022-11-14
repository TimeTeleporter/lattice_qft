#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    lattice::Lattice3d,
    simulation::{ComputationResult, ComputationType3d, Observable, Simulation3d, TestSim3d},
    REL_TEMP_ARY,
};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const _SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

const TEST_RANGE: usize = 8;

// Parameters for the Wilson loop
const _TEMP: f64 = 0.1;
const WIDTH: usize = 1;
const HEIGHT: usize = 1;
const BURNIN: usize = 10_000;
const ITERATIONS: usize = 100_000;

const TEST_PATH: &str = "data/test_results.csv";
const RESULTS_PATH: &str = "./data/sim_results.csv";
const BINNING_PATH: &str = "./data/sim_binning.csv";

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    // Initialise the simulations
    let mut comps: Vec<ComputationType3d<TEST_X, TEST_Y, TEST_T>> = Vec::new();
    for temp in REL_TEMP_ARY {
        let name: String =
            format!("{TEST_X}x{TEST_Y}x{TEST_T} Test for a {WIDTH}x{HEIGHT} Wilson loop");
        comps.push(ComputationType3d::Test(TestSim3d::new(
            name,
            Observable::Wilson(WIDTH, HEIGHT, temp),
            &lattice,
            temp,
            TEST_RANGE,
        )));
        let name: String =
            format!("{TEST_X}x{TEST_Y}x{TEST_T} Metropolis Sim for a {WIDTH}x{HEIGHT} Wilson loop");
        comps.push(ComputationType3d::Simulation(Simulation3d::new(
            name,
            lattice_qft::simulation::Algorithm::Metropolis,
            Observable::Wilson(WIDTH, HEIGHT, temp),
            &lattice,
            temp,
            BURNIN,
            ITERATIONS,
        )));
        let name: String =
            format!("{TEST_X}x{TEST_Y}x{TEST_T} Cluster Sim for a {WIDTH}x{HEIGHT} Wilson loop");
        comps.push(ComputationType3d::Simulation(Simulation3d::new(
            name,
            lattice_qft::simulation::Algorithm::Cluster,
            Observable::Wilson(WIDTH, HEIGHT, temp),
            &lattice,
            temp,
            BURNIN,
            ITERATIONS,
        )));
    }

    // Parallel over all temp data
    let data: Vec<ComputationResult> = comps
        .par_iter()
        .map(|comp| {
            let result = comp.run();
            result
        })
        .collect();

    for entry in data {
        if let Err(err) = entry.write_to_csv(RESULTS_PATH, BINNING_PATH, TEST_PATH) {
            eprint!("TestData Error: {}", err);
        };
    }
}
