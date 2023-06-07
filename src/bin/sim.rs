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
    outputdata::{return_correlation_lengths, Observe, OutputData, Plotting},
};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

const MAX_X: usize = 20;
const MAX_Y: usize = 20;
const MAX_T: usize = 20;
const SIZE: usize = MAX_X * MAX_Y * MAX_T;

const _RANGE: usize = 12;

const BURNIN: usize = 10_000;
const ITERATIONS: usize = 1_000_000;

/// The measurements for the wilson loop
const _WIDTH: usize = MAX_X / 3;
const _HEIGHT: usize = MAX_T;

const RESULTS_PATH: &str = "./data/results.csv";
const PLOT_PATH_INCOMPLETE: &str = "./data/plot_data/";

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    // Initialise the simulations
    let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
    for temp in lattice_qft::TEMP_ARY {
        let mut observables: Vec<OutputData<3, SIZE>> = Vec::new();
        observables.push(OutputData::new_action_observable(temp));
        //observables.push(OutputData::new_difference_plot(&lattice));
        //observables.push(OutputData::new_energy_plot(&lattice));
        observables.push(OutputData::new_correlation_plot(&lattice));
        comps.push(Computation::new_simulation(
            &lattice,
            temp,
            AlgorithmType::new_metropolis(),
            BURNIN,
            ITERATIONS,
            observables.clone(),
        ));
        /*
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
        comps.push(Computation::new_simulation(
            &lattice,
            temp,
            AlgorithmType::new_cluster(),
            BURNIN,
            ITERATIONS,
            observables.clone(),
        ));*/
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

        let (mut summary, outputs) = ComputationSummary::from_computation(computation, index);

        for output in outputs.into_iter() {
            match output {
                OutputData::ActionObservable(obs) => {
                    summary = summary.set_action(obs.result(), None);
                }
                OutputData::TestActionObservable(obs) => {
                    summary = summary.set_action(obs.result(), None);
                }
                OutputData::DifferencePlot(obs) => {
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
                    summary = summary.set_bonds_data();
                }
                OutputData::EnergyPlot(obs) => {
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
                    summary = summary.set_energy_data();
                }
                OutputData::CorrelationData(obs) => {
                    if let Some(plot) = obs.plot().into_iter().next() {
                        let (correlation_length12, correlation_length23, correlation_lenght13) =
                            return_correlation_lengths(plot.clone(), MAX_T);

                        let path: &str = &(PLOT_PATH_INCOMPLETE.to_owned()
                            + &"correlation_"
                            + &index.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path) {
                            eprint!("{}", err);
                        };

                        summary = summary.set_correlation_data().set_correlation_lenght(
                            correlation_length12,
                            correlation_length23,
                            correlation_lenght13,
                        );
                    }
                }
            }
        }

        if let Err(err) = summary.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        };
    }
}
