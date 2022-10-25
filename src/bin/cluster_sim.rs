#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]
#![recursion_limit = "2048"]

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    cluster::cluster_simulation3d,
    export::{calculate_binned_error, clean_csv, BinData, CsvData, SimResult},
    lattice::Lattice3d,
};

const MAX_X: usize = 100;
const MAX_Y: usize = 100;
const MAX_T: usize = 100;

const RESULTS_PATH: &str = "data/cluster_sim/cluster_results.csv";
const DATA_PATH: &str = "data/cluster_sim/cluster_data.csv";

const BURNIN: usize = 10_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 1_000_000;

const STACK_SIZE: usize = 1024 * 1024 * 1024;

fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(22)
        .stack_size(STACK_SIZE)
        .build_global()
        .unwrap();

    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let results: Vec<(SimResult, Vec<BinData>)> = [0.1]
        .par_iter()
        //.filter(|&&temp| temp == 0.001)
        .map(|&temp| {
            let (mut result, data) = cluster_simulation3d(&lattice, temp, BURNIN, ITERATIONS);
            let bins: Vec<BinData> = data.calculate_binnings(temp, 3);
            let error = match calculate_binned_error(&bins) {
                Ok(error) => Some(error),
                Err(err) => {
                    eprintln!("{err}");
                    None
                }
            };
            result.set_error(error);
            println!("{:?}", result);
            (result, bins)
        })
        .collect();

    for (res, bins) in results {
        if let Err(err) = res.read_write_csv(RESULTS_PATH) {
            eprint!("{err}");
        }

        if let Err(err) = clean_csv(DATA_PATH) {
            eprint!("{err}");
        };

        for bin in bins {
            if let Err(err) = bin.read_write_csv(DATA_PATH) {
                eprint!("{err}");
            }
        }
    }

    println!("Data stored");
}
