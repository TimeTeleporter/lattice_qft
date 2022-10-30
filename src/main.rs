#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use std::{ops::Deref, time::Instant};

use rayon::prelude::{IntoParallelIterator, ParallelIterator};

use lattice_qft::{
    error::BinDataAry,
    export::{clean_csv, CsvData},
    lattice::Lattice3d,
    simulation::{ClusterSim, MetrplsSim, SimResult, Simulation, Simulations3d},
    CLUSTER_BINNING_PATH, CLUSTER_RESULTS_PATH, METRPLS_BINNING_PATH, METRPLS_RESULTS_PATH,
};

const MAX_X: usize = 5;
const MAX_Y: usize = 5;
const MAX_T: usize = 5;

const BURNIN: usize = 10_000; // Number of sweeps until it starts counting.
const ITERATIONS: usize = 10_000_000;

fn main() {
    // Initialize the lattice
    let lattice: Lattice3d<MAX_X, MAX_Y, MAX_T> = Lattice3d::new();

    let temp: f64 = 0.1;

    let mut sims: Vec<(
        Simulations3d<MAX_X, MAX_Y, MAX_T, BURNIN, ITERATIONS>,
        (&str, &str),
    )> = Vec::new();

    let name: String = format!("{MAX_X}x{MAX_Y}x{MAX_T} ClusterSim");
    let paths: (&str, &str) = (CLUSTER_RESULTS_PATH, CLUSTER_BINNING_PATH);
    let sim: ClusterSim<3, { MAX_X * MAX_Y * MAX_T }, BURNIN, ITERATIONS> =
        ClusterSim::new(name, lattice.deref(), temp);
    sims.push((Simulations3d::ClusterSim3d(sim), paths));

    let name: String = format!("{MAX_X}x{MAX_Y}x{MAX_T} MetrplsSim");
    let paths: (&str, &str) = (METRPLS_RESULTS_PATH, METRPLS_BINNING_PATH);
    let metrpls_sim: MetrplsSim<3, { MAX_X * MAX_Y * MAX_T }, BURNIN, ITERATIONS> =
        MetrplsSim::new(name, lattice.deref(), temp);
    sims.push((Simulations3d::MetropolisSim3d(metrpls_sim), paths));

    let results: Vec<((SimResult, BinDataAry), (&str, &str))> = sims
        .into_par_iter()
        .map(|(sim, paths)| {
            #[allow(unused_mut)]
            let (mut result, data) = match sim {
                Simulations3d::ClusterSim3d(mut sim) => {
                    let time = Instant::now();
                    let result = sim.run();
                    let duration: f32 = time.elapsed().as_secs_f32();
                    println!("The Cluster Simulation took {duration} seconds.");
                    result
                }
                Simulations3d::MetropolisSim3d(mut sim) => {
                    let time = Instant::now();
                    let result = sim.run();
                    let duration: f32 = time.elapsed().as_secs_f32();
                    println!("The Metropolis Simulation took {duration} seconds.");
                    result
                }
            };
            let bins: BinDataAry = data.calculate_binnings(result.temp, 3);
            /*let error = bins
                .clone()
                .calculate_binned_error_tanh()
                .or_else(|err| {
                    eprintln!("{err}");
                    Err(err)
                })
                .ok();
            result.set_error(error);*/
            ((result, bins), paths)
        })
        .collect();

    for ((res, bins), (results_path, binning_path)) in results {
        if let Err(err) = res.read_write_csv(results_path) {
            eprint!("{err}");
        };
        if let Err(err) = clean_csv(binning_path) {
            eprint!("{err}");
        };
        for bin in bins {
            if let Err(err) = bin.read_write_csv(binning_path) {
                eprint!("{err}");
            }
        }
    }
}
