#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use lattice_qft::{
    error::BinDataAry,
    export::{clean_csv, CsvData},
    lattice::Lattice3d,
    simulation::{SimResult, Simulation3d, SimulationType},
};

const CUBE: usize = 2;

const MAX_X: usize = CUBE;
const MAX_Y: usize = CUBE;
const MAX_T: usize = CUBE;

const BURNIN: usize = 10_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 1_000_000;

const TEMP: f64 = 0.1;

const RESULTS_PATH: &str = "./data/sim_results.csv";
const BINNING_PATH: &str = "./data/sim_binning.csv";

fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let mut sims: Vec<Simulation3d<MAX_X, MAX_Y, MAX_T>> = Vec::new();
    let name: String = format!("{MAX_X}x{MAX_Y}x{MAX_T} Cluster Simulation, Size normalized");
    let size_normalized: bool = true;
    sims.push(Simulation3d::new(
        name,
        SimulationType::ClusterSim,
        size_normalized,
        &lattice,
        TEMP,
        BURNIN,
        ITERATIONS,
    ));
    let name: String = format!("{MAX_X}x{MAX_Y}x{MAX_T} Metropolis Simulation, Size normalized");
    sims.push(Simulation3d::new(
        name,
        SimulationType::MetropolisSim,
        size_normalized,
        &lattice,
        TEMP,
        BURNIN,
        ITERATIONS,
    ));

    let results: Vec<(SimResult, BinDataAry)> = sims
        .into_par_iter()
        .map(|sim| {
            let (result, data) = sim.run();
            let bins: BinDataAry = data.calculate_binnings(result.temp, 3);
            (result, bins)
        })
        .collect();

    for (res, bins) in results {
        if let Err(err) = res.read_write_csv(RESULTS_PATH) {
            eprint!("{err}");
        };
        if let Err(err) = clean_csv(BINNING_PATH) {
            eprint!("{err}");
        };
        for bin in bins {
            if let Err(err) = bin.read_write_csv(BINNING_PATH) {
                eprint!("{err}");
            }
        }
    }
}

/*
let error = bins.clone().calculate_binned_error_tanh().or_else(|err| {
    eprintln!("{err}");
    Err(err)
    }).ok();
    result.set_error(error);
*/
