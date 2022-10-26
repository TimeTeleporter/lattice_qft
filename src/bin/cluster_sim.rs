#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use std::time::Instant;

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    export::{calculate_binned_error, clean_csv, BinData, CsvData},
    lattice::Lattice3d,
    simulation::{cluster_simulation3d, SimResult},
    CLUSTER_BINNING_PATH, CLUSTER_RESULTS_PATH,
};

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

const BURNIN: usize = 10_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 1_000_000;

fn main() {
    let time = Instant::now();
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

    println!("Sim done, time elapsed: {} s", time.elapsed().as_secs());

    for (res, bins) in results {
        if let Err(err) = res.read_write_csv(CLUSTER_RESULTS_PATH) {
            eprint!("{err}");
        }

        if let Err(err) = clean_csv(CLUSTER_BINNING_PATH) {
            eprint!("{err}");
        };

        for bin in bins {
            if let Err(err) = bin.read_write_csv(CLUSTER_BINNING_PATH) {
                eprint!("{err}");
            }
        }
    }

    println!("Data stored, time elapsed: {} s", time.elapsed().as_secs());
}
