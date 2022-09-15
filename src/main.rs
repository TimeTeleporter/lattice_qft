#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
//#![allow(dead_code)]

use field3d::Field3d;
use lattice3d::Lattice3d;
use observable::Action;

#[allow(unused_imports)]
use metropolis::Metropolis;

#[allow(unused_imports)]
use cluster::Cluster;

mod cluster;
mod field3d;
mod lattice3d;
mod metropolis;
mod observable;

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    const TEST_X: usize = 2;
    const TEST_Y: usize = 2;
    const TEST_T: usize = 2;

    const TEST_RANGE: usize = 8;
    const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

    const PERMUTATIONS: usize = TEST_RANGE.pow(SIZE as u32); // 8 ^ 8 = 16_777_216

    const EQUIL: usize = 100_000; // Number of sweeps until it starts counting.
    const MAX_TRIES: usize = 1_000_000_000;
    const EPSILON: f64 = 0.001;

    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    // Saving all possible configurations
    let mut configurations: Vec<Field3d<i8, TEST_X, TEST_Y, TEST_T>> = Vec::new();
    let mut field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);
    configurations.push(field.clone());

    for _ in 0..PERMUTATIONS {
        'updateconfig: for index in 0..SIZE {
            match field.values[index] {
                0 | 1 | 2 | 3 | 4 | 5 | 6 => {
                    field.values[index] = field.values[index] + 1;
                    break 'updateconfig;
                }
                7 => {
                    field.values[index] = 0;
                }
                _ => {
                    panic!("config entry out of bounds.");
                }
            }
        }
        configurations.push(field.clone());
    }

    let mut partfn: f64 = 0.0;
    for field in configurations.into_iter() {
        partfn = partfn + (-field.action_observable()).exp();
    }
    let partfn = partfn / PERMUTATIONS as f64;
    println!("Calculated all configurations.");

    // Initialize a field to compare against
    let field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::random(&lattice);
    field.print_values_formated();
    let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field);

    // Sweeps to achieve equilibrium
    println!("Started equilibrium phase.");
    for _ in 0..EQUIL {
        field.metropolis_sweep();
    }

    let mut test: f64 = 0.0;
    for sweeps in 0..MAX_TRIES {
        field.metropolis_sweep();
        test = test + field.action_observable();
        if sweeps % 1000 == 0 {
            println!("Z: {}, <Z>: {}", partfn, test / sweeps as f64);
        }
        if ((test / sweeps as f64) - partfn).abs() < EPSILON {
            return;
        }
    }
    panic!("Didnt achieve equilibrium");
}
