use std::ops::{Add, Mul, Sub};

use rand::prelude::*;

use crate::field3d::Field3d;

const VERBOSE: bool = false;

pub trait Metropolis3d {
    type FieldType: Add<Output = Self::FieldType>
        + Default
        + Sub<Output = Self::FieldType>
        + From<i8>
        + Copy
        + Mul<Output = Self::FieldType>;

    /// Counts the number of lattice points, depends on the implementation;s
    const LATTICEPOINTS: usize;

    const TEMP: f64;

    fn metropolis_single(&mut self, index: usize, rng: &mut ThreadRng);

    fn metropolis_random(&mut self) {
        let mut rng = ThreadRng::default();
        let index: usize = rng.gen_range(0..Self::LATTICEPOINTS);
        self.metropolis_single(index, &mut rng);
    }

    fn metropolis_sweep(&mut self) {
        let mut rng = ThreadRng::default();
        for index in 0..Self::LATTICEPOINTS {
            self.metropolis_single(index, &mut rng);
        }
    }

    fn calculate_action(&self) -> i64;

    fn calculate_link_action(site: Self::FieldType, neighbour: Self::FieldType) -> Self::FieldType {
        site * site - site * neighbour * 2_i8.into() + neighbour * neighbour
    }
}

/// Central place to change the implemented type.
type Int = i32;

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Metropolis3d
    for Field3d<'a, Int, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type FieldType = Int;

    const LATTICEPOINTS: usize = MAX_X * MAX_Y * MAX_T;

    const TEMP: f64 = 0.01;

    fn metropolis_single(&mut self, index: usize, rng: &mut ThreadRng) {
        // Initialize the actions to be comapred
        let value = self.get_value(index).clone();
        let coin: bool = rng.gen();
        let new_value = match coin {
            true => value.clone() + 1_i8 as Self::FieldType,
            false => value.clone() - 1_i8 as Self::FieldType,
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
        let prob: f64 = (f64::from(action - new_action) * Self::TEMP).exp();
        if draw <= prob {
            self.values[index] = new_value;
            if VERBOSE {
                println!(
                    "Accepted: prob: {:.4}, draw: {:.4}, action: {}, new_action: {}",
                    prob, draw, action, new_action
                );
            }
        } else if VERBOSE {
            println!(
                "Declined: prob: {:.4}, draw: {:.4}, action: {}, new_action: {}",
                prob, draw, action, new_action
            );
        }
    }

    fn calculate_action(&self) -> i64 {
        let mut action: i64 = 0;
        for index in 0..Self::LATTICEPOINTS {
            let value = self.values[index];
            for neighbour in self.lattice.get_neighbours_array(index) {
                let neighbour = self.values[neighbour];
                action = action + Self::calculate_link_action(value, neighbour) as i64;
            }
        }
        action
    }
}

#[test]
fn test_action_add_one() {
    use crate::lattice3d::Lattice3d;

    const TEST_X: usize = 10;
    const TEST_Y: usize = 10;
    const TEST_T: usize = 10;

    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();
    let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);

    let action = field.calculate_action();

    for value in field.values.iter_mut() {
        *value = *value + 1;
    }

    assert_eq!(action, field.calculate_action());
}
