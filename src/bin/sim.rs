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
    outputdata::{ObservableOutputData, Observe, OutputData, PlotOutputData, Plotting},
};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

const TEST_X: usize = 50;
const TEST_Y: usize = 50;
const TEST_T: usize = 50;
const SIZE: usize = TEST_X * TEST_Y * TEST_T;

const _RANGE: usize = 12;

const BURNIN: usize = 10_000;
const ITERATIONS: usize = 100_000;

/// The measurements for the wilson loop
const WIDTH: usize = TEST_X / 5;
const HEIGHT: usize = TEST_T;

const RESULTS_PATH: &str = "./data/results.csv";
const PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/";

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
        observables.push(OutputData::new_difference_plot(&lattice));
        observables.push(OutputData::new_energy_plot(&lattice));
        observables.push(OutputData::new_correlation_plot(&lattice));
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
    }

    // Parallel over all temp data
    let data: Vec<Computation<3, SIZE>> = comps
        .into_par_iter()
        .filter_map(|comp| comp.run().ok())
        .collect();

    println!("All calcualtions done");

    for (_, computation) in data.into_iter().enumerate() {
        let index: usize =
            match ComputationSummary::fetch_csv_data(RESULTS_PATH).and_then(|summary| {
                summary
                    .last()
                    .map(|last| last.index + 1)
                    .ok_or("No last element, starting anew.".into())
            }) {
                Ok(index) => index,
                Err(err) => {
                    eprint!("{}", err);
                    0
                }
            };

        let (mut new, outputs) = ComputationSummary::from_computation(computation, index);

        for output in outputs.into_iter() {
            match output {
                OutputData::Observable(ObservableOutputData::Action(obs)) => {
                    new = new.set_action(obs.result(), None);
                }
                OutputData::Observable(ObservableOutputData::TestAction(obs)) => {
                    new = new.set_action(obs.result(), None);
                }
                OutputData::Plot(PlotOutputData::Difference(obs)) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(PLOT_PATH_INCOMPLETE.to_owned()
                            + &"difference_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path) {
                            eprint!("{}", err);
                        };
                    }
                    new = new.set_bonds_data();
                }
                OutputData::Plot(PlotOutputData::Energy(obs)) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(PLOT_PATH_INCOMPLETE.to_owned()
                            + &"energy_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path) {
                            eprint!("{}", err);
                        };
                    }
                    new = new.set_energy_data();
                }
                OutputData::Plot(PlotOutputData::Correlation(obs)) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(PLOT_PATH_INCOMPLETE.to_owned()
                            + &"correlation_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path) {
                            eprint!("{}", err);
                        };
                    }
                    new = new.set_correlation_data();
                }
            }
        }

        if let Err(err) = new.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        };
    }
}
