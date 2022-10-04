use serde::{Deserialize, Serialize};

use lattice_qft::{
    action::Action, export::CsvData, field::Field, lattice::Lattice, metropolis::Metropolis,
};

const DIM: usize = 1;
const SIZE_ARY: [usize; DIM] = [2];
const SIZE: usize = SIZE_ARY[0];

const TEMP: f64 = 1.0;
const BURNIN: usize = 100_000;
const ITERATIONS: usize = 1_000_000_000;

const PATH: &str = "data/line_sim/line_sim.csv";

/// Struct to print out the results of the line simulation.
#[derive(Debug, Serialize, Deserialize)]
struct LineResults {
    temp: f64,
    coupling: f64,
    burnin: usize,
    iterations: usize,
    observable: f64,
    variance: f64,
}

fn main() {
    let lattice: Lattice<DIM, SIZE> = Lattice::new(SIZE_ARY);

    println!("{:?}", lattice);

    let result: LineResults = line_sim(&lattice, TEMP);

    if let Err(err) = result.read_write_csv(PATH) {
        eprint!("Output Error: {}", err);
    }
}

fn line_sim(lattice: &Lattice<DIM, SIZE>, temp: f64) -> LineResults {
    let field: Field<i8, DIM, SIZE> = Field::new(lattice);
    let mut field: Field<i32, DIM, SIZE> = Field::from_field(field);

    for _ in 0..BURNIN {
        field.metropolis_sweep(temp);
    }

    let mut obs_ary: Vec<f64> = Vec::with_capacity(ITERATIONS);
    for _ in 0..ITERATIONS {
        field.metropolis_sweep(temp);
        obs_ary.push(field.lattice_action(temp));
    }

    let observable: f64 = obs_ary.iter().sum::<f64>() / ITERATIONS as f64;
    let variance: f64 = obs_ary
        .into_iter()
        .map(|x| (x - observable) * (x - observable))
        .sum::<f64>()
        / ITERATIONS as f64;

    let coupling: f64 = (TEMP * 2.0).sqrt();

    LineResults {
        coupling,
        temp,
        burnin: BURNIN,
        iterations: ITERATIONS,
        observable,
        variance,
    }
}
