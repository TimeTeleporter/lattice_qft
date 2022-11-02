use std::ops::{Add, Deref, DerefMut, Div, Mul, Sub};

use rand::prelude::*;

use crate::field::{Field, Field3d};

pub trait Action {
    type FieldType: Copy
        + Add<Output = Self::FieldType>
        + Sub<Output = Self::FieldType>
        + Mul<Output = Self::FieldType>
        + Div<Output = Self::FieldType>
        + Default
        + From<i8>;

    const D: usize;
    const SIZE: usize;

    /// Sums up the action across all lattice bonds
    fn sum_link_actions(&self) -> i64;

    /// Calculates the action of the lattice with the coupling constant.
    fn action_observable(&self, temp: f64) -> f64 {
        (self.sum_link_actions() as f64) * temp
    }

    /// Calculates the [action_observable], normalized by the amount of lattice sites.
    fn size_normalized_action_observable(&self, temp: f64) -> f64 {
        self.action_observable(temp) / Self::SIZE as f64
    }

    /// Calculates the [action_observable], normalized by the amount of lattice sites.
    fn bond_normalized_action_observable(&self, temp: f64) -> f64 {
        self.action_observable(temp) / (Self::D * Self::SIZE) as f64
    }

    /// Calculates ```(x - y)^2``` for two values.
    fn calculate_link_action(site: Self::FieldType, neighbour: Self::FieldType) -> Self::FieldType {
        (site - neighbour) * (site - neighbour)
    }

    fn calculate_assumed_action(&self, index: usize, site: Self::FieldType) -> i64;

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
    const D: usize = D;
    const SIZE: usize = SIZE;

    fn sum_link_actions(&self) -> i64 {
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

    fn calculate_assumed_action(&self, index: usize, value: Self::FieldType) -> i64 {
        let mut action: i64 = 0;
        for neighbour in self.lattice.get_neighbours_array(index) {
            let neighbour = self.get_value(neighbour).clone();
            action = action + Self::calculate_link_action(value, neighbour).into();
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
    const D: usize = 3;
    const SIZE: usize = MAX_X * MAX_Y * MAX_T;

    /// Sums up all link actions of the lattice.
    fn sum_link_actions(&self) -> i64 {
        self.deref().sum_link_actions()
    }

    fn calculate_assumed_action(&self, index: usize, site: Self::FieldType) -> i64 {
        self.deref().calculate_assumed_action(index, site)
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
    use crate::lattice::Lattice3d;

    let lattice: Lattice3d<9, 9, 9> = Lattice3d::new();

    let field: Field3d<i8, 9, 9, 9> = Field3d::random(&lattice);

    let field: Field3d<i32, 9, 9, 9> = Field3d::from_field(field);

    assert!(field.sum_link_actions() >= 0);
}

#[test]
fn test_action_new_zero() {
    use crate::lattice::Lattice;

    let lattice: Lattice<4, 81> = Lattice::new([3, 3, 3, 3]);
    let field: Field<i8, 4, 81> = Field::new(&lattice);

    assert_eq!(field.sum_link_actions(), 0);
}

#[test]
fn test_action_given_field() {
    use crate::lattice::Lattice3d;

    let lattice: Lattice3d<3, 3, 3> = Lattice3d::new();
    let mut field: Field3d<i8, 3, 3, 3> = Field3d::new(&lattice);

    field.values[0] = 1;

    assert_eq!(field.sum_link_actions(), 6);

    field.values[4] = 1;

    assert_eq!(field.sum_link_actions(), 12);

    field.values[1] = 2;

    assert_eq!(field.sum_link_actions(), 28)
}

#[test]
fn test_action_add_one() {
    use crate::lattice::Lattice3d;

    const TEST_X: usize = 10;
    const TEST_Y: usize = 10;
    const TEST_T: usize = 10;

    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::new();
    let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);

    let action = field.sum_link_actions();

    for value in field.values.iter_mut() {
        *value = *value + 1;
    }

    assert_eq!(action, field.sum_link_actions());
}

#[test]
fn test_random_normalize() {
    use crate::lattice::Lattice3d;

    let lattice: Lattice3d<2, 2, 2> = Lattice3d::new();
    let field: Field3d<i8, 2, 2, 2> = Field3d::random(&lattice);
    let mut field: Field3d<i32, 2, 2, 2> = Field3d::from_field(field);

    field.normalize_random();

    assert!(field.values.iter().find(|&&x| x == 0).is_some());
}

#[test]
fn test_calculate_assumed_action() {
    use crate::lattice::Lattice3d;

    let lattice: Lattice3d<10, 10, 10> = Lattice3d::new();
    let field: Field3d<i8, 10, 10, 10> = Field3d::new(&lattice);
    let mut field: Field3d<i32, 10, 10, 10> = Field3d::from_field(field);

    let index: usize = 234;
    let value: i32 = 3;

    let assumed: i64 = field.calculate_assumed_action(index, value);

    field.values[index] = 3;

    let test: i64 = field.sum_link_actions();

    println!("test: {test}, assumed: {assumed}");

    assert_eq!(test, assumed);
}
