use std::ops::{Add, DerefMut, Div, Mul, Sub};

use rand::prelude::*;

use crate::{
    action::Action,
    field::{Field, Field3d},
};

pub trait Metropolis: Action {
    fn metropolis_single(&mut self, index: usize, temp: f64, rng: &mut ThreadRng);

    fn metropolis_random(&mut self, temp: f64) {
        let mut rng = ThreadRng::default();
        let index: usize = rng.gen_range(0..Self::SIZE);
        self.metropolis_single(index, temp, &mut rng);
    }

    fn metropolis_sweep(&mut self, temp: f64) {
        let mut rng = ThreadRng::default();
        for index in 0..Self::SIZE {
            self.metropolis_single(index, temp, &mut rng);
        }
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Metropolis for Field<'a, T, D, SIZE>
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
    [(); D * 2_usize]:,
{
    fn metropolis_single(&mut self, index: usize, temp: f64, rng: &mut ThreadRng) {
        // Initialize the change to be measured
        let mut new_field = self.clone();
        let coin: bool = rng.gen();
        new_field.values[index] = match coin {
            true => new_field.values[index] + Self::FieldType::from(1_i8),
            false => new_field.values[index] - Self::FieldType::from(1_i8),
        };

        // Calculate the action of both possibilities
        let action = self.action_observable();
        let new_action = new_field.action_observable();

        // Accept the new action if its lower than the previous.
        // Else accept it with a proportional probability.
        let draw: f64 = rng.gen_range(0.0..=1.0);
        let prob: f64 = (((action - new_action) as f64) * temp).exp();
        if draw <= prob {
            *self = new_field;
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
    fn metropolis_single(&mut self, index: usize, temp: f64, rng: &mut ThreadRng) {
        self.deref_mut().metropolis_single(index, temp, rng)
    }
}
