use rand::{distributions::Standard, prelude::Distribution, random};

use crate::lattice3d::{Directions, Lattice3d};

pub struct Field3d<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub values: Vec<T>,
    pub lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
    //updated: Vec<usize>,
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'_, T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn get_prev_neighbour_value(&self, index: usize, direction: Directions) -> &T {
        self.get_value(self.lattice.get_prev_neighbour_index(index, direction))
    }

    pub fn get_next_neighbour_value(&self, index: usize, direction: Directions) -> &T {
        self.get_value(self.lattice.get_next_neighbour_index(index, direction))
    }

    pub fn get_value_from_coordinates(&self, coordinates: (usize, usize, usize)) -> &T {
        self.get_value(self.lattice.get_index_from_coordinates(coordinates))
    }

    pub fn get_value(&self, index: usize) -> &T {
        &self.values[index]
    }

    #[cfg(test)]
    #[allow(dead_code)]
    pub fn set_value_from_coordinates(&mut self, value: T, coordinates: (usize, usize, usize)) {
        self.set_value(value, self.lattice.get_index_from_coordinates(coordinates));
    }

    pub fn set_value(&mut self, value: T, index: usize) {
        self.values[index] = value;
        //self.updated.push(index);
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'_, T, MAX_X, MAX_Y, MAX_T>
where
    T: std::fmt::Debug,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Prints the value in a intuitive way.
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
    Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    #[allow(dead_code)]
    /// Constructor for a new field on the lattice initialized to be zero everywhere.
    pub fn new(lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>) -> Self {
        let values: Vec<()> = vec![(); MAX_X * MAX_Y * MAX_T];
        let values: Vec<T> = values.iter().map(|_| T::default()).collect();

        Field3d::<'a, T, MAX_X, MAX_Y, MAX_T> {
            values,
            lattice,
            //updated: Vec::<usize>::new(),
        }
    }
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    Standard: Distribution<T>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Constructor for a new field on the lattice initialized to be random everywhere.
    pub fn random(lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>) -> Self {
        let values: Vec<()> = vec![(); MAX_X * MAX_Y * MAX_T];
        let values: Vec<T> = values.iter().map(|_| random()).collect();

        Field3d::<'a, T, MAX_X, MAX_Y, MAX_T> {
            values,
            lattice,
            //updated: Vec::<usize>::new(),
        }
    }
}

pub trait ConvertField<'a, T, U, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    U: From<T>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn from_field(field: Field3d<'a, T, MAX_X, MAX_Y, MAX_T>) -> Self;
}

impl<'a, T, U, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    ConvertField<'a, T, U, MAX_X, MAX_Y, MAX_T> for Field3d<'a, U, MAX_X, MAX_Y, MAX_T>
where
    U: From<T>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn from_field(field: Field3d<'a, T, MAX_X, MAX_Y, MAX_T>) -> Self {
        let values = field.values.into_iter().map(|x| U::from(x)).collect();

        Field3d::<'a, U, MAX_X, MAX_Y, MAX_T> {
            values,
            lattice: field.lattice,
        }
    }
}

/// Tests the constructor, setting a value and printing formatted.
#[test]
fn test_field() {
    const TEST_X: usize = 2;
    const TEST_Y: usize = 2;
    const TEST_T: usize = 2;
    type TestInt = i32;

    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::default();

    let mut field: Field3d<TestInt, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);

    let value = 3;

    field.set_value_from_coordinates(value, (2, 2, 1));

    assert_eq!(field.values, vec![0, 0, 0, 0, value, 0, 0, 0]);
}

#[test]
fn test_field_conversion() {
    const TEST_X: usize = 2;
    const TEST_Y: usize = 2;
    const TEST_T: usize = 2;

    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::default();

    let field8: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);

    let field16: Field3d<i16, TEST_X, TEST_Y, TEST_T> = Field3d::new(&lattice);

    let field: Field3d<i16, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field8);

    assert_eq!(field16.values, field.values);
}
