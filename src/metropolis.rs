use std::ops::{Add, Div, Mul, Sub};

use crate::{field3d::Field3d, observable::Action};
use rand::prelude::*;

pub trait Metropolis: Action {
    fn metropolis_single(&mut self, index: usize, rng: &mut ThreadRng);

    fn metropolis_random(&mut self) {
        let mut rng = ThreadRng::default();
        let index: usize = rng.gen_range(0..Self::SIZE);
        self.metropolis_single(index, &mut rng);
    }

    fn metropolis_sweep(&mut self) {
        let mut rng = ThreadRng::default();
        for index in 0..Self::SIZE {
            self.metropolis_single(index, &mut rng);
        }
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Metropolis
    for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    f64: From<T>,
    T: Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Default
        + From<i8>
        + Into<i64>
        + Into<f64>
        + PartialOrd
        + Copy,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn metropolis_single(&mut self, index: usize, rng: &mut ThreadRng) {
        // Initialize the actions to be comapred
        let value = self.get_value(index).clone();
        let coin: bool = rng.gen();
        let new_value = match coin {
            true => value.clone() + Self::FieldType::from(1_i8),
            false => value.clone() - Self::FieldType::from(1_i8),
        };

        // Calculate the actions
        let mut action = Self::FieldType::default();
        let mut new_action = Self::FieldType::default();
        for neighbour in self.lattice.get_neighbours_array(index) {
            let neighbour = self.get_value(neighbour).clone();
            action = action + Self::calculate_link_action(value, neighbour);
            new_action = new_action + Self::calculate_link_action(new_value, neighbour);
        }

        // Accept the new action if its lower than the previous.
        // Else accept it with a proportional probability.
        let draw: f64 = rng.gen_range(0.0..1.0);
        let prob: f64 = (Self::TEMP * f64::from(action - new_action)).exp();
        if draw <= prob {
            self.values[index] = new_value;
        }
    }
}

/// We initialize a 2 by 2 by 2 lattice, on which all possible configurations
/// with values from 0 to 8 are known. Then we run a metropolis simulation
/// in order to test that it converges to the desired distribution.
#[test]
fn test_all_configurations() {
    use crate::lattice3d::Lattice3d;

    const TEST_X: usize = 2;
    const TEST_Y: usize = 2;
    const TEST_T: usize = 2;

    const TEST_RANGE: usize = 8;
    const SIZE: usize = TEST_X * TEST_Y * TEST_Y; // 8 lattice points

    const PERMUTATIONS: usize = TEST_RANGE.pow(SIZE as u32); // 8 ^ 8 = 16_777_216

    const MAXTRIES: usize = 1_000_000_000;
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

    // Initialize a field to compare against
    let field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::random(&lattice);
    let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field);

    let mut test: f64 = 0.0;
    for sweeps in 0..MAXTRIES {
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
