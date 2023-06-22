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
    outputdata::{Observe, OutputData, Plotting},
};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

const CUBE: usize = 54;

const MAX_X: usize = CUBE;
const MAX_Y: usize = CUBE;
const MAX_T: usize = CUBE;
const SIZE: usize = MAX_X * MAX_Y * MAX_T;

const _RANGE: usize = CUBE;

const BURNIN: usize = 100_000;
const ITERATIONS: usize = 1_000_000;

/// The measurements for the wilson loop
const _WIDTH: usize = MAX_X / 3;
const _HEIGHT: usize = MAX_T;

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    // Initialise the simulations
    let mut comps: Vec<Computation<3, SIZE>> = Vec::new();
    for temp in lattice_qft::INVESTIGATE_ARY3_LOWER {
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
    }

    // Parallel over all temp data
    let data: Vec<Computation<3, SIZE>> = comps
        .into_par_iter()
        .filter_map(|comp| comp.run().ok())
        .collect();

    println!("All calcualtions done");

    // Parse the results
    for (_, computation) in data.into_iter().enumerate() {
        // Fetching the index to append to
        let index: usize = match ComputationSummary::fetch_csv_data(lattice_qft::RESULTS_PATH, true)
            .and_then(|summary| {
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

        // Building the computation summary and handling outputs.
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
                        let path: &str = &(lattice_qft::PLOT_PATH_INCOMPLETE.to_owned()
                            + &"difference_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path, true) {
                            eprint!("{}", err);
                        };
                    }
                    summary = summary.set_bonds_data();
                }
                OutputData::EnergyPlot(obs) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(lattice_qft::PLOT_PATH_INCOMPLETE.to_owned()
                            + &"energy_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path, true) {
                            eprint!("{}", err);
                        };
                    }
                    summary = summary.set_energy_data();
                }
                OutputData::CorrelationData(obs) => {
                    if let Some(plot) = obs.plot().into_iter().next() {
                        let path: &str = &(lattice_qft::PLOT_PATH_INCOMPLETE.to_owned()
                            + &"correlation_"
                            + &index.to_string()
                            + &".csv");
                        if let Err(err) = plot.overwrite_csv(path) {
                            eprint!("{}", err);
                        };
                        summary = summary.set_correlation_data();
                    }
                }
            }
        }

        if let Err(err) = summary.read_write_csv(lattice_qft::RESULTS_PATH, true) {
            eprint!("{}", err);
        };
    }
}
