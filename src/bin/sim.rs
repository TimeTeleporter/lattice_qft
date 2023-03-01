#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use lattice_qft::{
    //REL_TEMP_ARY,
    algorithm::AlgorithmType,
    computation::{Computation, ComputationResult, Compute},
    export::{clean_csv, CsvData},
    lattice::{Lattice, Lattice3d},
    observable::ObservableType,
    plot::PlotType,
};

const TEST_X: usize = 10;
const TEST_Y: usize = 10;
const TEST_T: usize = 10;
const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points
const PLOTSIZE: usize = TEST_X * TEST_Y;

const _RANGE: usize = 16;

const BURNIN: usize = 10_000;
const ITERATIONS: usize = 1_000_000;

const WIDTH: usize = 1;
const HEIGHT: usize = 1;

const RESULTS_PATH: &str = "./data/results.csv";
const PLOT_PATH: &str = "./data/plot.csv";

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
        comps.push(Computation::new_wilson_sim(
            &lattice,
            temp,
            AlgorithmType::new_metropolis(),
            BURNIN,
            ITERATIONS,
            WIDTH,
            HEIGHT,
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
        let (export, opt) = entry.into_export();

        if let Some(plot) = opt {
            if let Err(err) = clean_csv(PLOT_PATH) {
                eprint!("{}", err);
            };

            for point in plot.values {
                if let Err(err) = point.read_write_csv(PLOT_PATH) {
                    eprint!("{}", err);
                };
            }
        }

        if let Err(err) = export.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        };
    }
}
