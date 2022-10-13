#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    cluster::cluster_simulation3d,
    export::{CsvData, SimResult},
    lattice::Lattice3d,
    LONG_TEMP_ARY,
};

const MAX_X: usize = 2;
const MAX_Y: usize = 2;
const MAX_T: usize = 2;

const RESULTS_PATH: &str = "data/cluster_sim/cluster_data.csv";

const BURNIN: usize = 100_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 10_000_000;

fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let results: Vec<SimResult> = LONG_TEMP_ARY
        .par_iter()
        //.filter(|&&temp| temp == 0.001)
        .map(|&temp| {
            let (result, _) = cluster_simulation3d(&lattice, temp, BURNIN, ITERATIONS);
            result
        })
        .collect();

    for res in results {
        if let Err(err) = res.read_write_csv(RESULTS_PATH) {
            eprint!("{}", err);
        }
    }

    println!("Done");
}
