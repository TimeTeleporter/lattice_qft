#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use lattice_qft::{
    computation::{
        Algorithm, Computatable, Computation, ComputationResults, Observable, Simulation,
    },
    export::CsvData,
    lattice::Lattice3d,
};

const CUBE: usize = 2;

const MAX_X: usize = CUBE;
const MAX_Y: usize = CUBE;
const MAX_T: usize = CUBE;
const SIZE: usize = MAX_X * MAX_Y * MAX_T;

const WIDTH: usize = 1;
const HEIGHT: usize = 1;

const BURNIN: usize = 100_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 10_000_000;

const TEMP: f64 = 0.1;

const RESULTS_PATH: &str = "./data/sim_results.csv";
//const BINNING_PATH: &str = "./data/sim_binning.csv";

fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let observable: Observable = Observable::Wilson {
        width: WIDTH,
        height: HEIGHT,
    };

    let mut sims: Vec<Computation<3, SIZE>> = Vec::new();
    sims.push(Simulation::new_compuatation(
        &lattice,
        Algorithm::Cluster,
        observable,
        TEMP,
        BURNIN,
        ITERATIONS,
    ));
    sims.push(Simulation::new_compuatation(
        &lattice,
        Algorithm::Metropolis,
        observable,
        TEMP,
        BURNIN,
        ITERATIONS,
    ));

    let results: Vec<ComputationResults> = sims
        .into_par_iter()
        .filter_map(|sim| sim.run().ok())
        .collect();

    for res in results {
        if let Err(err) = res.read_write_csv(RESULTS_PATH) {
            eprint!("{err}");
        };
    }
}

/*
let error = bins.clone().calculate_binned_error_tanh().or_else(|err| {
    eprintln!("{err}");
    Err(err)
    }).ok();
    result.set_error(error);
*/
