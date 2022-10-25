use std::ops::{Add, Deref, DerefMut, Mul, Sub};

use rand::{distributions::Standard, prelude::Distribution, random};

use crate::lattice::{Lattice, Lattice3d, LatticeCoords};

#[derive(Debug, Clone)]
pub struct Field<'a, T, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    pub values: Vec<T>,
    pub lattice: &'a Lattice<D, SIZE>,
}

// Implementing getting values
impl<T, const D: usize, const SIZE: usize> Field<'_, T, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn get_value_from_coords(&self, coords: LatticeCoords<D>) -> &T {
        self.get_value(self.lattice.calc_index_from_coords(coords))
    }

    pub fn get_value(&self, index: usize) -> &T {
        &self.values[index]
    }
}

// Implementing shift
impl<T, const D: usize, const SIZE: usize> Field<'_, T, D, SIZE>
where
    T: Copy + Sub<Output = T>,
    [(); D * 2_usize]:,
{
    /// Subtracts a shift from all values of the field.
    pub fn shift_values(&mut self, shift: T) {
        for value in self.values.iter_mut() {
            *value = *value - shift;
        }
    }
}

// Implementing mirroring
impl<T, const D: usize, const SIZE: usize> Field<'_, T, D, SIZE>
where
    T: Copy + From<i8> + Add<Output = T> + Sub<Output = T> + Mul<Output = T>,
    [(); D * 2_usize]:,
{
    /// Mirrors all field values on a plane with a modifier
    pub fn mirror_values(mut self, plane: T, modifier: T) -> Self {
        for value in self.values.iter_mut() {
            *value = <i8 as Into<T>>::into(2_i8) * plane + modifier - *value;
        }
        self
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Field<'a, T, D, SIZE>
where
    T: Default,
    [(); D * 2_usize]:,
{
    #[allow(dead_code)]
    /// Constructor for a new field on the lattice initialized to be zero everywhere.
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        let values: Vec<()> = vec![(); lattice.size.iter().product()];
        let values: Vec<T> = values.iter().map(|_| T::default()).collect();

        Field::<'a, T, D, SIZE> { values, lattice }
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Field<'a, T, D, SIZE>
where
    T: Default + PartialOrd + Copy,
    [(); D * 2_usize]:,
{
    pub fn get_max(&self) -> (usize, T) {
        let mut result: (usize, T) = (usize::default(), T::default());
        for (index, &item) in self.values.iter().enumerate() {
            if item > result.1 {
                result = (index, item);
            }
        }
        result
    }

    pub fn get_min(&self) -> (usize, T) {
        let mut result: (usize, T) = (usize::default(), T::default());
        for (index, &item) in self.values.iter().enumerate() {
            if item < result.1 {
                result = (index, item);
            }
        }
        result
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Field<'a, T, D, SIZE>
where
    Standard: Distribution<T>,
    [(); D * 2_usize]:,
{
    #[allow(dead_code)]
    /// Constructor for a new field on the lattice initialized to be random everywhere.
    pub fn random(lattice: &'a Lattice<D, SIZE>) -> Self {
        let values: Vec<()> = vec![(); lattice.size.iter().product()];
        let values: Vec<T> = values.iter().map(|_| random()).collect();

        Field::<'a, T, D, SIZE> { values, lattice }
    }
}

impl<'a, U, const D: usize, const SIZE: usize> Field<'a, U, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Convert a field of one type into a field of a different type losslessly.
    pub fn from_field<T: Into<U>>(field: Field<'a, T, D, SIZE>) -> Self {
        let values = field.values.into_iter().map(|x| x.into()).collect();

        Field::<'a, U, D, SIZE> {
            values,
            lattice: field.lattice,
        }
    }
}

impl<T, const SIZE: usize> Field<'_, T, 3, SIZE>
where
    T: std::fmt::Debug,
{
    /// Prints the value in a intuitive way.
    pub fn print_values_formated(&self, size: [usize; 3]) {
        for t in 0..size[2] {
            println!("t = {}", t);
            for y in 0..size[1] {
                print!("[");
                for x in 0..size[0] {
                    if x == size[0] - 1 {
                        println!(
                            "{:?} ]",
                            self.get_value_from_coords(LatticeCoords::new([x, y, t]))
                        );
                    } else {
                        print!(
                            "{:?}, ",
                            self.get_value_from_coords(LatticeCoords::new([x, y, t]))
                        );
                    }
                }
            }
        }
    }
}

impl<T, const SIZE: usize> Field<'_, T, 2, SIZE>
where
    T: std::fmt::Debug,
{
    /// Prints the value in a intuitive way.
    pub fn print_values_formated(&self, size: [usize; 2]) {
        for y in 0..size[1] {
            print!("[");
            for x in 0..size[0] {
                if x == size[0] - 1 {
                    println!(
                        "{:?} ]",
                        self.get_value_from_coords(LatticeCoords::new([x, y]))
                    );
                } else {
                    print!(
                        "{:?}, ",
                        self.get_value_from_coords(LatticeCoords::new([x, y]))
                    );
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------

/// Implements a filed of given Type T on a 3-dimensional lattice of size
/// MAX_X * MAX_Y * MAX_T.
#[derive(Debug, Clone)]
pub struct Field3d<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    Field<'a, T, 3, { MAX_X * MAX_Y * MAX_T }>,
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

// By implementing deref we can utilize all the functions of the tuple element
// without implementing it ourselves.
impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Deref
    for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type Target = Field<'a, T, 3, { MAX_X * MAX_Y * MAX_T }>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

// Here we also need to implement DerefMut, as we will need to change values of the set field.
impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> DerefMut
    for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'_, T, MAX_X, MAX_Y, MAX_T>
where
    T: std::fmt::Debug,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Prints the value in a intuitive way
    pub fn print_values_formated(&self) {
        self.0.print_values_formated([MAX_X, MAX_Y, MAX_T]);
    }
}

impl<'a, T: 'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    #[allow(dead_code)]
    /// Constructor for a new field on the lattice initialized to be zero everywhere.
    pub fn new(lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>) -> Self {
        Field3d::<'a, T, MAX_X, MAX_Y, MAX_T>(Field::<'a, T, 3, { MAX_X * MAX_Y * MAX_T }>::new(
            lattice,
        ))
    }
}

impl<'a, T: 'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    Standard: Distribution<T>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    #[allow(dead_code)]
    /// Constructor for a new field on the lattice initialized to be random everywhere.
    pub fn random(lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>) -> Self {
        Field3d::<'a, T, MAX_X, MAX_Y, MAX_T>(Field::<T, 3, { MAX_X * MAX_Y * MAX_T }>::random(
            lattice,
        ))
    }
}

impl<'a, U, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Field3d<'a, U, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Convert a field of one type into a field of a different type losslessly.
    pub fn from_field<T: Into<U>>(field: Field3d<'a, T, MAX_X, MAX_Y, MAX_T>) -> Self {
        Field3d::<'a, U, MAX_X, MAX_Y, MAX_T>(Field::<U, 3, { MAX_X * MAX_Y * MAX_T }>::from_field(
            field.0,
        ))
    }
}

#[test]
fn test_field_conversion() {
    let lattice: Lattice<3, 8> = Lattice::new([2, 2, 2]);

    let field8: Field<i8, 3, 8> = Field::new(&lattice);

    let field16: Field<i16, 3, 8> = Field::new(&lattice);

    let field: Field<i16, 3, 8> = Field::from_field(field8);

    assert_eq!(field16.values, field.values);
}

#[test]
fn test_field_shift() {
    let lattice: Lattice<4, 81> = Lattice::new([3, 3, 3, 3]);

    let mut field: Field<i8, 4, 81> = Field::new(&lattice);

    field.shift_values(3);

    assert_eq!(field.values[80], -3);
}

#[test]
fn test_field_mirror() {
    use crate::action::Action;
    let lattice: Lattice<4, 81> = Lattice::new([3, 3, 3, 3]);
    let field: Field<i8, 4, 81> = Field::random(&lattice);
    let mut field: Field<i32, 4, 81> = Field::from_field(field);

    let old_action: i64 = field.action_observable();

    field = field.mirror_values(29, 0);

    assert_eq!(old_action, field.action_observable());
}
