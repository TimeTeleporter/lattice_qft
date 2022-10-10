use std::ops::{Add, Div, Mul, Sub};

use crate::{action::Action, field::Field};

use rand::{rngs::ThreadRng, Rng};

pub trait Cluster: Action {
    fn cluster_sweep(&mut self, temp: f64);
}

impl<'a, T, const D: usize, const SIZE: usize> Cluster for Field<'a, T, D, SIZE>
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
    fn cluster_sweep(&mut self, temp: f64) {
        let mut rng = ThreadRng::default();

        // Set the mirror plane randomly on a height value
        let plane: T = self.values[rng.gen_range(0..SIZE)];
        let modifier: T = match rng.gen::<bool>() {
            true => 0_i8.into(),
            fales => match rng.gen::<bool>() {
                true => (-1_i8).into(),
                false => 1_i8.into(),
            },
        };

        let mut reflected: Field<T, D, SIZE> = self.clone();
        reflected.mirror_values(plane, modifier);
        let reflected: Field<T, D, SIZE> = reflected;
    }
}
