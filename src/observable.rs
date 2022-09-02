use std::ops::{Add, Mul, Sub};

use crate::field3d::{Field, Field3d};

pub trait Action {
    type FieldType: Copy
        + Add<Output = Self::FieldType>
        + Sub<Output = Self::FieldType>
        + Mul<Output = Self::FieldType>
        + Default
        + From<i8>;

    const SIZE: usize;

    fn calculate_action(&self) -> i64;

    fn calculate_link_action(site: Self::FieldType, neighbour: Self::FieldType) -> Self::FieldType {
        site * site - site * neighbour * 2_i8.into() + neighbour * neighbour
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Action for Field<'a, T, D, SIZE>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Default + From<i8> + Into<i64>,
    [(); D * 2_usize]:,
{
    type FieldType = T;
    const SIZE: usize = SIZE;

    fn calculate_action(&self) -> i64 {
        let mut action: i64 = 0;
        for index in 0..Self::SIZE {
            let value = self.values[index];
            for neighbour in self.lattice.get_neighbours_array(index) {
                let neighbour = self.values[neighbour];
                action = action + Self::calculate_link_action(value, neighbour).into();
            }
        }
        action / Self::SIZE as i64
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Action
    for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Default + From<i8> + Into<i64>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type FieldType = T;
    const SIZE: usize = MAX_X * MAX_Y * MAX_T;

    fn calculate_action(&self) -> i64 {
        let mut action: i64 = 0;
        for index in 0..Self::SIZE {
            let value = self.values[index];
            for neighbour in self.lattice.get_neighbours_array(index) {
                let neighbour = self.values[neighbour];
                action = action + Self::calculate_link_action(value, neighbour).into();
            }
        }
        action / Self::SIZE as i64
    }
}
