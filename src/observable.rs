use std::ops::{Add, Div, Mul, Sub};

use rand::prelude::*;

use crate::field3d::{Field, Field3d};

pub trait Action {
    type FieldType: Copy
        + Add<Output = Self::FieldType>
        + Sub<Output = Self::FieldType>
        + Mul<Output = Self::FieldType>
        + Div<Output = Self::FieldType>
        + Default
        + From<i8>;

    const SIZE: usize;

    fn calculate_action(&self) -> i64;

    fn calculate_link_action(site: Self::FieldType, neighbour: Self::FieldType) -> Self::FieldType {
        site * site - site * neighbour * 2_i8.into() + neighbour * neighbour
    }

    fn normalize(&mut self);
    fn normalize_random(&mut self);
}

impl<'a, T, const D: usize, const SIZE: usize> Action for Field<'a, T, D, SIZE>
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Default
        + From<i8>
        + Into<i64>
        + PartialOrd,
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

    fn normalize(&mut self) {
        let shift = (self.get_max().1 + self.get_min().1) / 2_i8.into();
        self.shift_values(shift);
    }

    fn normalize_random(&mut self) {
        let mut rng = ThreadRng::default();
        let &shift = self.values.choose(&mut rng).unwrap();
        self.shift_values(shift);
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Action
    for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Default
        + From<i8>
        + Into<i64>
        + PartialOrd,
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

    fn normalize(&mut self) {
        let shift = (self.get_max().1 + self.get_min().1) / 2_i8.into();
        self.shift_values(shift);
    }

    fn normalize_random(&mut self) {
        let mut rng = ThreadRng::default();
        let &shift = self.values.choose(&mut rng).unwrap();
        self.shift_values(shift);
    }
}

#[test]
fn test_action_non_negative() {
    use crate::lattice3d::Lattice3d;

    let lattice: Lattice3d<9, 9, 9> = Lattice3d::new();

    let field: Field3d<i8, 9, 9, 9> = Field3d::random(&lattice);

    let field: Field3d<i32, 9, 9, 9> = Field3d::from_field(field);

    assert!(field.calculate_action() >= 0);
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

#[test]
fn test_random_normalize() {
    use crate::lattice3d::Lattice3d;

    let lattice: Lattice3d<2, 2, 2> = Lattice3d::new();
    let field: Field3d<i8, 2, 2, 2> = Field3d::random(&lattice);
    let mut field: Field3d<i32, 2, 2, 2> = Field3d::from_field(field);

    field.normalize_random();

    assert_eq!(field.values.iter().find(|&&x| x == 0).is_some(), true);
}
