use crate::{field3d::Field3d, observable::Action};
use rand::prelude::*;

pub trait Metropolis: Action {
    const TEMP: f64;

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

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Metropolis
    for Field3d<'a, i32, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
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
        }
    }
}
