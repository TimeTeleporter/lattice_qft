#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]
#![feature(split_array)]

use lattice_qft::{
    export::CsvData, field3d::Field3d, lattice3d::Lattice3d, metropolis::Metropolis,
    observable::Action,
};

use num::complex::ComplexFloat;
use serde::{Deserialize, Serialize};

const TEST_X: usize = 2;
const TEST_Y: usize = 2;
const TEST_T: usize = 2;
const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

const TEST_PATH: &str = "data/test_data16.csv";
const SIM_PATH: &str = "data/sim_data16.csv";

/// Datatype to save and read simulation output.
#[derive(Debug, Serialize, Deserialize)]
struct TestData {
    temp: f64,
    test_range: usize,
    test_observable: f64,
}

/// Datatype to save and read simulation output.
#[derive(Debug, Serialize, Deserialize)]
struct SimData {
    temp: f64,
    sim_burnin: usize,
    sim_iterations: usize,
    sim_observable: f64,
}

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
fn main() {
    const EXPO_ARY: [f64; 41] = [
        -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6,
        -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    ];

    // Initialize the lattice
    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();

    for expo in EXPO_ARY {
        let temp = expo.expf(10.0);

        if let Err(err) = TestData::calculate_configurations(&lattice, temp).write_to_csv(TEST_PATH)
        {
            eprint!("TestData Error: {}", err);
        };

        println!("{}: TestData calculated", temp);

        if let Err(err) = SimData::metropolis_simulation(&lattice, temp).write_to_csv(SIM_PATH) {
            eprint!("SimData Error: {}", err);
        };

        println!("{}: Done", temp);
    }
}

impl SimData {
    fn metropolis_simulation(lattice: &Lattice3d<TEST_X, TEST_Y, TEST_T>, temp: f64) -> SimData {
        const BURNIN: usize = 10_000; // Number of sweeps until it starts counting.
        const MAX_TRIES: usize = 1_000_000;

        // Initialize a field to compare against
        let field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::random(&lattice);
        //field.print_values_formated();
        let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field);

        // Sweeps to achieve equilibrium
        for _ in 0..BURNIN {
            field.metropolis_sweep(temp);
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        let mut sim_observable: f64 = 0.0;
        for _sweeps in 0..MAX_TRIES {
            field.metropolis_sweep(temp);
            sim_observable = sim_observable + field.lattice_action(temp);
        }
        let sim_observable: f64 = sim_observable / MAX_TRIES as f64;

        SimData {
            temp,
            sim_burnin: BURNIN,
            sim_iterations: MAX_TRIES,
            sim_observable,
        }
    }
}

impl TestData {
    fn calculate_configurations(
        lattice: &Lattice3d<TEST_X, TEST_Y, TEST_T>,
        temp: f64,
    ) -> TestData {
        const TEST_RANGE: usize = 16;
        const PERMUTATIONS: usize = TEST_RANGE.pow(SIZE as u32 - 1); // 16 ^ 7 = 268’435’456
        const BOUDARY: i8 = TEST_RANGE as i8 - 1;

        let mut field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::new(lattice);

        field.values[7] = (TEST_RANGE / 2) as i8;

        let mut test_ary: Vec<f64> = Vec::with_capacity(PERMUTATIONS);
        let mut bolz_ary: Vec<f64> = Vec::with_capacity(PERMUTATIONS);

        for _ in 0..PERMUTATIONS {
            'updateconfig: for index in 0..SIZE {
                match field.values[index] {
                    x if x < BOUDARY => {
                        field.values[index] = field.values[index] + 1;
                        break 'updateconfig;
                    }
                    x if x == BOUDARY => {
                        field.values[index] = 0;
                    }
                    _ => {
                        panic!("config entry out of bounds.");
                    }
                }
            }
            let field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field.clone());
            //field.print_values_formated();
            let bolz: f64 = (-field.lattice_action(temp)).exp();
            test_ary.push(field.lattice_action(temp) * bolz);
            bolz_ary.push(bolz);
        }

        // The partition function is the sum over all Bolzmann weights
        let partfn: f64 = bolz_ary.into_iter().sum();

        // The observable is the sum over all weights (Bolzmann times observable),
        // devided by the partition function.
        let test: f64 = test_ary.into_iter().sum();
        TestData {
            temp,
            test_range: TEST_RANGE,
            test_observable: test / partfn,
        }
    }
}
