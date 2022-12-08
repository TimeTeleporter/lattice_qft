//! This module implements a methods on fields that are of a HeightVariable type.

use std::{
    fmt::Display,
    ops::{Add, Mul, Sub},
};

use rand::prelude::*;

use crate::field::Field;

/// A collection of traits to be fullfilled by the field type, i.e. i8, i32
pub trait HeightVariable<T>:
    Copy
    + Add<Output = T>
    + Sub<Output = T>
    + Mul<Output = T>
    + Default
    + From<i8>
    + Into<f64>
    + Display
    + PartialOrd
{
}

impl<T> HeightVariable<T> for T where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Default
        + From<i8>
        + Into<f64>
        + Display
        + PartialOrd
{
}

/// For [Field]s of [HeightVariable]s we implement possible methods.
pub trait HeightField<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    /// Calculates the difference of two height variables.
    fn difference(&self, site1: usize, site2: usize) -> T;

    /// Calculates the square difference of two height variables.
    fn difference_squared(&self, site1: usize, site2: usize) -> T {
        let diff: T = self.difference(site1, site2);
        diff * diff
    }

    /// Sums up all bond differences connected to a site.
    fn local_bond_diff_square_sum(&self, index: usize) -> T;

    /// Sums up all bond differences in positive coordinate direction
    /// connected to a site.
    fn pos_bond_diff_square_sum(&self, index: usize) -> T;

    /// Sums up all bond differences.
    fn global_bond_diff_square_sum(&self) -> T {
        let mut sum: T = T::default();
        for index in 0..SIZE {
            sum = sum + self.pos_bond_diff_square_sum(index);
        }
        sum
    }

    /// Shifts all values by the value of a random chosen height value.
    fn normalize_random(&mut self);

    /// Subtracts a shift from all values of the field.
    fn shift_values(&mut self, shift: T);

    fn reflect_value(value: T, plane: T, modifier: T) -> T {
        Into::<T>::into(2_i8) * plane + modifier - value
    }
}

impl<'a, T, const D: usize, const SIZE: usize> HeightField<T, D, SIZE> for Field<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn difference(&self, site1: usize, site2: usize) -> T {
        self.values[site1] - self.values[site2]
    }

    fn local_bond_diff_square_sum(&self, index: usize) -> T {
        let mut sum: T = T::default();
        for neighbour in self.lattice.get_neighbours_array(index) {
            sum = sum + self.difference_squared(index, neighbour);
        }
        sum
    }

    fn pos_bond_diff_square_sum(&self, index: usize) -> T {
        let mut sum: T = T::default();
        for neighbour in self.lattice.pos_neighbours_array(index) {
            sum = sum + self.difference_squared(index, neighbour);
        }
        sum
    }

    fn normalize_random(&mut self) {
        let mut rng = ThreadRng::default();
        let &shift = self.values.choose(&mut rng).unwrap();
        self.shift_values(shift);
    }

    /// Subtracts a shift from all values of the field.
    fn shift_values(&mut self, shift: T) {
        for value in self.values.iter_mut() {
            *value = *value - shift;
        }
    }
}
