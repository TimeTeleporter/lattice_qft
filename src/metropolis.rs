use crate::{field3d::Field3d, observable::Action};
use rand::prelude::*;
use std::ops::{Add, Div, Mul, Sub};

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
        // Initialize the change to be measured
        let mut new_field = self.clone();
        let coin: bool = rng.gen();
        new_field.values[index] = match coin {
            true => new_field.values[index] + Self::FieldType::from(1_i8),
            false => new_field.values[index] - Self::FieldType::from(1_i8),
        };

        // Calculate the action of both possibilities
        let action = self.lattice_action(temp);
        let new_action = new_field.lattice_action(temp);

        // Accept the new action if its lower than the previous.
        // Else accept it with a proportional probability.
        let draw: f64 = rng.gen_range(0.0..1.0);
        let prob: f64 = (action - new_action).exp();
        if draw <= prob {
            *self = new_field;
        }
    }
}

/*
        // Initialize the actions to be comapred
        let value = self.get_value(index).clone();
        let coin: bool = rng.gen();
        let new_value = match coin {
            true => value.clone() + Self::FieldType::from(1_i8),
            false => value.clone() - Self::FieldType::from(1_i8),
        };

        // Calculate the action of both possibilities
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
*/
