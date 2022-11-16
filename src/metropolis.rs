use std::ops::{Add, DerefMut, Div, Mul, Sub};

use rand::prelude::*;

use crate::{
    action::Action,
    field::{Field, Field3d},
};

pub trait Metropolis<const D: usize, const SIZE: usize>: Action<D, SIZE> {
    fn metropolis_single(&mut self, index: usize, temp: f64, rng: &mut ThreadRng) -> bool;

    fn metropolis_random(&mut self, temp: f64) {
        let mut rng = ThreadRng::default();
        let index: usize = rng.gen_range(0..SIZE);
        self.metropolis_single(index, temp, &mut rng);
    }

    /// Performs a metropolis step for each lattice site. Returns the number
    /// of accepted metropolis steps.
    fn metropolis_sweep(&mut self, temp: f64) -> usize {
        let mut rng = ThreadRng::default();
        let mut acceptance: usize = 0;
        for index in 0..SIZE {
            if self.metropolis_single(index, temp, &mut rng) {
                acceptance += 1;
            };
        }
        acceptance
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Metropolis<D, SIZE> for Field<'a, T, D, SIZE>
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
        + Copy
        + Sync,
    [(); D * 2_usize]:,
{
    fn metropolis_single(&mut self, index: usize, temp: f64, rng: &mut ThreadRng) -> bool {
        // Initialize the change to be measured
        let coin: bool = rng.gen();
        let new_value = match coin {
            true => self.values[index] + Self::FieldType::from(1_i8),
            false => self.values[index] - Self::FieldType::from(1_i8),
        };

        // Calculate the action of both possibilities
        let action = self.calculate_assumed_action(index, self.values[index]);
        let new_action = self.calculate_assumed_action(index, new_value);

        // Accept the new action if its lower than the previous.
        // Else accept it with a proportional probability.
        let draw: f64 = rng.gen_range(0.0..=1.0);
        let prob: f64 = (((action - new_action) as f64) * temp).exp();
        if draw <= prob {
            self.values[index] = new_value;
            return true;
        }
        false
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Metropolis<3, { MAX_X * MAX_Y * MAX_T }> for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
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
        + Copy
        + Sync,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn metropolis_single(&mut self, index: usize, temp: f64, rng: &mut ThreadRng) -> bool {
        self.deref_mut().metropolis_single(index, temp, rng)
    }
}
