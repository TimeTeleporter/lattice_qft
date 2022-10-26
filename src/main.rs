#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use std::time::Instant;

use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use lattice_qft::{
    error::{BinDataAry, ObsChain},
    export::{clean_csv, CsvData},
    lattice::Lattice3d,
    simulation::{cluster_simulation3d, metropolis_simulation3d, SimResult},
    CLUSTER_BINNING_PATH, CLUSTER_RESULTS_PATH, METRPLS_BINNING_PATH, METRPLS_RESULTS_PATH,
};

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

const BURNIN: usize = 10_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 1_000_000;

fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let sims: [(
        &str,
        &str,
        fn(&Lattice3d<MAX_X, MAX_Y, MAX_T>, f64, usize, usize) -> (SimResult, ObsChain),
    ); 2] = [
        (
            CLUSTER_RESULTS_PATH,
            CLUSTER_BINNING_PATH,
            metropolis_simulation3d,
        ),
        (
            METRPLS_RESULTS_PATH,
            METRPLS_BINNING_PATH,
            cluster_simulation3d,
        ),
    ];

    let _res: Vec<()> = sims
        .par_iter()
        .map(|(result_path, binning_path, sim_fn)| {
            process(&lattice, result_path, binning_path, *sim_fn)
        })
        .collect();
}

fn process(
    lattice: &Lattice3d<MAX_X, MAX_Y, MAX_T>,
    result_path: &str,
    binning_path: &str,
    simulation: fn(&Lattice3d<MAX_X, MAX_Y, MAX_T>, f64, usize, usize) -> (SimResult, ObsChain),
) {
    let time = Instant::now();

    let results: Vec<(SimResult, BinDataAry)> = [0.1]
        .par_iter()
        .map(|&temp| {
            let (mut result, data) = simulation(&lattice, temp, BURNIN, ITERATIONS);
            let bins: BinDataAry = data.calculate_binnings(temp, 3);
            let error = bins
                .clone()
                .calculate_binned_error_tanh()
                .or_else(|err| {
                    eprintln!("{err}");
                    Err(err)
                })
                .ok();
            result.set_error(error);
            println!("{:?}", result);
            (result, bins)
        })
        .collect();

    println!("Sim done, time elapsed: {} s", time.elapsed().as_secs());

    for (res, bins) in results {
        if let Err(err) = res.read_write_csv(result_path) {
            eprint!("{err}");
        }

        if let Err(err) = clean_csv(binning_path) {
            eprint!("{err}");
        };

        for bin in bins {
            if let Err(err) = bin.read_write_csv(binning_path) {
                eprint!("{err}");
            }
        }
    }

    println!("Data stored, time elapsed: {} s", time.elapsed().as_secs());
}
