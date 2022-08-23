use std::ops::Sub;

use crate::lattice3d::{Directions, Lattice3d};

pub struct ZeroForm<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    values: Vec<T>,
    lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
    updated: Vec<usize>,
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    ZeroForm<'_, T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn get_value_from_coordinates(&self, coordinates: (usize, usize, usize)) -> &T {
        self.get_value(self.lattice.get_index_from_coordinates(coordinates))
    }

    pub fn get_value(&self, index: usize) -> &T {
        &self.values[index]
    }

    pub fn _set_value_from_coordinates(&mut self, value: T, coordinates: (usize, usize, usize)) {
        self.set_value(value, self.lattice.get_index_from_coordinates(coordinates));
    }

    pub fn set_value(&mut self, value: T, index: usize) {
        self.values[index] = value;
        self.updated.push(index);
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    ZeroForm<'_, T, MAX_X, MAX_Y, MAX_T>
where
    T: Copy,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn set_line_to_value(
        &mut self,
        value: T,
        direction: Directions,
        start_coordinates: (usize, usize, usize),
        end_coordinates: (usize, usize, usize),
    ) {
        let start_index: usize = self.lattice.get_index_from_coordinates(start_coordinates);
        let end_index: usize = self.lattice.get_index_from_coordinates(end_coordinates);
        let mut index = start_index;

        loop {
            self.set_value(value, index);
            index = self.lattice.get_next_neighbour_index(index, direction);
            if index == start_index || index == end_index {
                self.set_value(value, index);
                return;
            }
        }
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    ZeroForm<'_, T, MAX_X, MAX_Y, MAX_T>
where
    T: std::fmt::Debug,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn print_values_formated(&self) {
        for t in 0..MAX_T {
            println!("t = {}", t);
            for y in 0..MAX_Y {
                print!("[");
                for x in 0..MAX_X {
                    if x == MAX_X - 1 {
                        println!("{:?} ]", self.get_value_from_coordinates((x, y, t)));
                    } else {
                        print!("{:?}, ", self.get_value_from_coordinates((x, y, t)));
                    }
                }
            }
        }
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    From<&'a Lattice3d<MAX_X, MAX_Y, MAX_T>> for ZeroForm<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn from(lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>) -> Self {
        let values = vec![(); MAX_X * MAX_Y * MAX_T];

        let values = values.iter().map(|_| T::default()).collect();

        ZeroForm::<'a, T, MAX_X, MAX_Y, MAX_T> {
            values,
            lattice,
            updated: Vec::<usize>::new(),
        }
    }
}

#[derive(Debug)]
/// Datatype to save the outgoing links in an array.
pub struct Link<T>([T; 3]);

impl<T> Default for Link<T>
where
    T: Default,
{
    fn default() -> Self {
        Link([(); 3].map(|_| T::default()))
    }
}

/// Datatype to set one 3-array per lattice site.
pub type OneForm<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> =
    ZeroForm<'a, Link<T>, MAX_X, MAX_Y, MAX_T>;

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    OneForm<'a, T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Overwrites all entries who or whose coderivatives have changed.
    /// This empties all 'updated' stacks.
    pub fn differential_update_all(
        &mut self,
        zero_form: &mut ZeroForm<'a, T, MAX_X, MAX_Y, MAX_T>,
    ) {
        let mut updated = Vec::<usize>::new();
        while let Some(index) = self.updated.pop() {
            updated.push(index);
        }
        while let Some(index) = zero_form.updated.pop() {
            updated.push(index);
        }
        for index in updated.into_iter() {
            self.differential_update_one(index);
        }
    }

    pub fn differential_update_one(&mut self, index: usize) {}
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    From<&ZeroForm<'a, T, MAX_X, MAX_Y, MAX_T>> for OneForm<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Default + Copy + Sub<Output = T>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Implements a transformation from a 0-form to a 1-form via the differential.
    fn from(zero_form: &ZeroForm<'a, T, MAX_X, MAX_Y, MAX_T>) -> Self {
        let mut one_form = OneForm::<T, MAX_X, MAX_Y, MAX_T>::from(zero_form.lattice);

        for iter in one_form.values.iter_mut().enumerate() {
            let (index, links): (usize, &mut Link<T>) = iter;
            for (i, link) in links.0.iter_mut().enumerate() {
                *link = zero_form.values[zero_form
                    .lattice
                    .get_next_neighbour_index(index, Directions::from(i))]
                    - zero_form.values[index];
            }
        }

        one_form
    }
}

/// Tests the constructor, setting a value and printing formatted.
#[test]
fn test_zero_form() {
    const TEST_X: usize = 2;
    const TEST_Y: usize = 2;
    const TEST_T: usize = 2;
    type TestInt = i32;

    let lattice = Lattice3d::<TEST_X, TEST_Y, TEST_T>::default();

    let mut zero_form = ZeroForm::<TestInt, TEST_X, TEST_Y, TEST_T>::from(&lattice);

    let value = 3;

    zero_form._set_value_from_coordinates(value, (2, 2, 1));

    assert_eq!(zero_form.values, vec![0, 0, 0, 0, value, 0, 0, 0]);
}

#[test]
fn test_zero_form_line() {
    const TEST_X: usize = 3;
    const TEST_Y: usize = 3;
    const TEST_T: usize = 3;
    type TestInt = i32;

    let lattice = Lattice3d::<TEST_X, TEST_Y, TEST_T>::default();

    let mut zero_form = ZeroForm::<TestInt, TEST_X, TEST_Y, TEST_T>::from(&lattice);

    let value = 3;

    zero_form.set_line_to_value(value, Directions::T, (0, 1, 2), (0, 1, 0));

    zero_form.print_values_formated();

    assert_eq!(
        zero_form.values,
        vec![
            0, 0, 0, value, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, value, 0, 0, 0, 0, 0
        ]
    );
}

#[test]
fn test_one_form_from_zero_form() {
    const TEST_X: usize = 100;
    const TEST_Y: usize = 100;
    const TEST_T: usize = 100;
    type TestInt = i32;

    let lattice = Lattice3d::<TEST_X, TEST_Y, TEST_T>::default();

    let zero_form = ZeroForm::<TestInt, TEST_X, TEST_Y, TEST_T>::from(&lattice);

    let mut one_form = OneForm::<TestInt, TEST_X, TEST_Y, TEST_T>::from(&zero_form);

    one_form._set_value_from_coordinates(Link([10, 10, 10]), (50, 50, 99))
}
