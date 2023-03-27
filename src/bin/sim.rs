#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]
#![feature(let_chains)]

use lattice_qft::{
    algorithm::AlgorithmType,
    computation::{Computation, ComputationSummary, Compute},
    export::CsvData,
    lattice::Lattice3d,
    outputdata::{ObservableOutputData, Observe, OutputData, PlotOutputData},
};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const SIZE: usize = TEST_X * TEST_Y * TEST_T;

const RANGE: usize = 16;

const BURNIN: usize = 10_000;
const ITERATIONS: usize = 100_000;

/// The measurements for the wilson loop
const WIDTH: usize = TEST_X / 2;
const HEIGHT: usize = TEST_T;

const RESULTS_PATH: &str = "./data/results.csv";
const _PLOT_PATH: &str = "./data/plot.csv";

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    // Initialise the simulations
    let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
    for temp in [0.1] {
        let mut observables: Vec<OutputData<3, SIZE>> = Vec::new();
        observables.push(OutputData::new_action_observable(temp));
        let mut test_observables: Vec<OutputData<3, SIZE>> = Vec::new();
        test_observables.push(OutputData::new_test_action_observable(temp));
        comps.push(Computation::new_simulation(
            &lattice,
            temp,
            AlgorithmType::new_metropolis(),
            BURNIN,
            ITERATIONS,
            observables.clone(),
        ));
        comps.push(Computation::new_wilson_sim(
            &lattice,
            temp,
            AlgorithmType::new_metropolis(),
            BURNIN,
            ITERATIONS,
            WIDTH,
            HEIGHT,
            observables,
        ));
        comps.push(Computation::new_test(
            &lattice,
            temp,
            RANGE,
            test_observables.clone(),
        ));
        comps.push(Computation::new_wilson_test(
            &lattice,
            temp,
            RANGE,
            WIDTH,
            HEIGHT,
            test_observables,
        ))
    }

    // Parallel over all temp data
    let data: Vec<Computation<3, SIZE>> = comps
        .into_par_iter()
        .filter_map(|comp| comp.run().ok())
        .collect();

    println!("All calcualtions done");

    for computation in data {
        let index: usize = match ComputationSummary::fetch_csv_data(RESULTS_PATH) {
            Ok(data) => {
                let mut index: usize = 0;
                if let Some(comp) = data.last() {
                    index = comp.index
                }
                index
            }
            Err(err) => {
                eprintln!("{err}");
                0
            }
        };

        let (mut new, outputs) = ComputationSummary::from_computation(computation, index);

        for output in outputs.into_iter() {
            new = match output {
                OutputData::Observable(ObservableOutputData::Action(obs)) => {
                    new.set_action(obs.result(), None)
                }
                OutputData::Observable(ObservableOutputData::TestAction(obs)) => {
                    new.set_action(obs.result(), None)
                }
                OutputData::Plot(PlotOutputData::Difference(_)) => todo!(),
                OutputData::Plot(PlotOutputData::Energy(_)) => todo!(),
            };
        }

        if let Err(err) = new.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        };
    }
}
