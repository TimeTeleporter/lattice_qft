use std::ops::{Add, Deref, DerefMut, Div, Mul, Sub};

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
    const TEMP: f64 = 400.0;

    fn action_observable(&self) -> i64;

    /// Calculates the action of the lattice with the coupling constant.
    fn lattice_action(&self) -> f64 {
        (self.action_observable() as f64) * Self::TEMP
    }

    /// Calculates ```(x - y)^2``` for two values.
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

    fn action_observable(&self) -> i64 {
        let mut action: i64 = 0;
        for index in 0..Self::SIZE {
            let value = self.get_value(index).clone();
            for neighbour in self.lattice.pos_neighbours_array(index) {
                let neighbour = self.get_value(neighbour).clone();
                action = action + Self::calculate_link_action(value, neighbour).into();
            }
        }
        action
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

    /// Sums up all link actions of the lattice.
    fn action_observable(&self) -> i64 {
        self.deref().action_observable()
    }

    fn normalize(&mut self) {
        self.deref_mut().normalize()
    }

    fn normalize_random(&mut self) {
        self.deref_mut().normalize_random()
    }
}

#[test]
fn test_action_non_negative() {
    use crate::lattice3d::Lattice3d;

    let lattice: Lattice3d<9, 9, 9> = Lattice3d::new();

    let field: Field3d<i8, 9, 9, 9> = Field3d::random(&lattice);

    let field: Field3d<i32, 9, 9, 9> = Field3d::from_field(field);

    assert!(field.action_observable() >= 0);
}

#[test]
fn test_action_new_zero() {
    use crate::lattice3d::Lattice;

    let lattice: Lattice<4, 81> = Lattice::new([3, 3, 3, 3]);
    let field: Field<i8, 4, 81> = Field::new(&lattice);

    assert_eq!(field.action_observable(), 0);
}

#[test]
fn test_action_given_field() {
    use crate::lattice3d::Lattice3d;

    let lattice: Lattice3d<3, 3, 3> = Lattice3d::new();
    let mut field: Field3d<i8, 3, 3, 3> = Field3d::new(&lattice);

    field.values[0] = 1;

    assert_eq!(field.action_observable(), 6);

    field.values[4] = 1;

    assert_eq!(field.action_observable(), 12);

    field.values[1] = 2;

    assert_eq!(field.action_observable(), 28)
}

#[test]
fn test_action_add_one() {
    use crate::lattice3d::Lattice3d;

    const TEST_X: usize = 10;
    const TEST_Y: usize = 10;
    const TEST_T: usize = 10;

    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();
    let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);

    let action = field.action_observable();

    for value in field.values.iter_mut() {
        *value = *value + 1;
    }

    assert_eq!(action, field.action_observable());
}

#[test]
fn test_random_normalize() {
    use crate::lattice3d::Lattice3d;

    let lattice: Lattice3d<2, 2, 2> = Lattice3d::new();
    let field: Field3d<i8, 2, 2, 2> = Field3d::random(&lattice);
    let mut field: Field3d<i32, 2, 2, 2> = Field3d::from_field(field);

    field.normalize_random();

    assert!(field.values.iter().find(|&&x| x == 0).is_some());
}
