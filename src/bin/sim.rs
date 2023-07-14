#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]
#![feature(let_chains)]

use lattice_qft::{
    algorithm::AlgorithmType,
    computation::{parse_simulation_results, Computation, Compute},
    lattice::Lattice3d,
    outputdata::OutputData,
};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

// Simulation parameters
const REPETITIONS: usize = 5;
const BURNIN: usize = 200_000;
const ITERATIONS: usize = 200_000;

// Lattice sizes (16, 24, 36, 54)
const CUBE: usize = 16;
const MAX_X: usize = CUBE;
const MAX_Y: usize = CUBE;
const MAX_T: usize = CUBE;
const SIZE: usize = MAX_X * MAX_Y * MAX_T;

// The measurements for the wilson loop
const _RANGE: usize = CUBE;
const _WIDTH: usize = MAX_X / 3;
const _HEIGHT: usize = MAX_T;

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Setting the global thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(16)
        .build_global()
        .unwrap();

    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    for rep in 1..=REPETITIONS {
        // Initialise the simulations
        let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
        println!("Starting rep {}", rep);
        for temp in [0.26; 15] {
            let mut observables: Vec<OutputData<3, SIZE>> = Vec::new();
            observables.push(OutputData::new_action_observable(temp));
            observables.push(
                OutputData::new_correlation_plot(&lattice)
                    .set_frequency(10)
                    .set_repetitions(10),
            );
            comps.push(Computation::new_simulation(
                &lattice,
                temp,
                AlgorithmType::new_cluster(),
                BURNIN,
                ITERATIONS,
                observables.clone(),
            ));
        }
        // Parallel over all temp data
        let data: Vec<Computation<3, SIZE>> = comps
            .clone()
            .into_par_iter()
            .filter_map(|comp| comp.run().ok())
            .collect();

        parse_simulation_results(data);
    }
}
