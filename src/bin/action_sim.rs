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
const REPETITIONS: u64 = 1;
const BURNIN: u64 = 0;
const ITERATIONS: u64 = 2_000;

// Machine parameters
const THREADS: usize = 21;

#[allow(unused_imports)]
use lattice_qft::INVESTIGATE_ARY10_54 as TEMP_ARY;

// Lattice sizes (16, 24, 36, 54)
const CUBE: usize = 54;
const MAX_X: usize = CUBE;
const MAX_Y: usize = CUBE;
const MAX_T: usize = CUBE;
const SIZE: usize = MAX_X * MAX_Y * MAX_T;

// The measurements for the wilson loop
const HEIGHT: usize = MAX_T;

fn main() {
    // Setting the global thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(THREADS)
        .build_global()
        .unwrap();

    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    for rep in 1..=REPETITIONS {
        // Initialise the simulations
        let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
        println!("Starting rep {}", rep);

        // Define the computations
        for temp in [0.2, 0.3] {
            let mut observables: Vec<OutputData<3, SIZE>> = Vec::new();
            observables.push(OutputData::new_action_observable(temp));

            comps.push(Computation::new_wilson_sim(
                &lattice,
                temp,
                AlgorithmType::new_metropolis(),
                BURNIN,
                ITERATIONS,
                MAX_X * 4 / 6,
                HEIGHT,
                observables.clone(),
            ));

            comps.push(Computation::new_simulation(
                &lattice,
                temp,
                AlgorithmType::new_metropolis(),
                BURNIN,
                ITERATIONS,
                observables.clone(),
            ));
        }

        // Parallel over all temp data
        let data: Vec<Computation<3, SIZE>> = comps
            .into_par_iter()
            .filter_map(|comp| comp.run().ok())
            .collect();

        parse_simulation_results(data);
    }
}
