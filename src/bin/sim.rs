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
    lattice::{Lattice, Lattice3d},
    observable::ObservableType,
    plot::PlotType,
};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points
const PLOTSIZE: usize = TEST_X * TEST_Y;

const RANGE: usize = 24;

const BURNIN: usize = 10_000;
const ITERATIONS: usize = 100_000_000;

const WIDTH: usize = 1;
const HEIGHT: usize = 1;

const RESULTS_PATH: &str = "./data/results.csv";

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();
    let plot_lattice: Lattice<2, PLOTSIZE> = Lattice::new([TEST_X, TEST_Y]);

    // Initialise the simulations
    let mut comps: Vec<Computation<3, SIZE, PLOTSIZE>> = Vec::new();
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
            Some(PlotType::ElectricField(&plot_lattice)),
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
            Some(PlotType::ElectricField(&plot_lattice)),
        ));
    }

    // Parallel over all temp data
    let data: Vec<ComputationResult<3, SIZE, PLOTSIZE>> = comps
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
