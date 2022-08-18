//! A three dimensional spacetime lattice saved to an array.
//!
//! The x coordinate is periodic from 0 to MAX_X, the same for y and t.
//!
//! Now consider the example of a <4,5,3> 3-dimensional lattice. We index the entries as follows:
//!
//!           y            t=0      y            t=1      y            t=2
//!         4 ^  16 17 18 19        ^  36 37 38 39        ^  56 57 58 59  
//!         3 |  12 13 14 15        |  32 33 34 35        |  52 53 54 55  
//!         2 |  8  9  10 11        |  28 29 30 31        |  48 49 50 51  
//!         1 |  4  5  6  7         |  24 25 26 27        |  44 45 46 47  
//!         0 |  0  1  2  3         |  20 21 22 23        |  40 41 42 43  
//!           .-----------> x       .-----------> x       .-----------> x
//!              0  1  2  3            0  1  2  3            0  1  2  3
//!
//! This example gets tested in the tests below.

/// As many directions as there are dimensions.
pub enum Direction {
    X,
    Y,
    T,
}

// This is the unoptimized lattice datatype. We optimize it further below.
#[derive(Debug)]
pub struct Lattice3d<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    pub [T; MAX_X * MAX_Y * MAX_T],
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Default
    for Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: Default,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn default() -> Lattice3d<T, MAX_X, MAX_Y, MAX_T> {
        Lattice3d([(); MAX_X * MAX_Y * MAX_T].map(|_| T::default()))
    }
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Lattice<MAX_X, MAX_Y, MAX_T>
    for Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
}

impl<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    LatticeValues<T, MAX_X, MAX_Y, MAX_T> for Lattice3d<T, MAX_X, MAX_Y, MAX_T>
where
    T: std::fmt::Debug,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    fn print_values(&self) {
        println!("{:?}", self.0);
    }

    fn get_value(&self, index: usize) -> &T {
        &self.0[index]
    }
}

/// Lattice operation trait concerning the conversion between indices and coordinates, implemented for three dimensions
pub trait Lattice<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> {
    fn prev_neighbour_index(&self, index: usize, direction: Direction) -> usize {
        let (mut x, mut y, mut t) = Self::get_coordinates_from_index(index);
        match direction {
            Direction::X => x = (x + MAX_X - 1) % MAX_X,
            Direction::Y => y = (y + MAX_Y - 1) % MAX_Y,
            Direction::T => t = (t + MAX_T - 1) % MAX_T,
        }

        Self::get_index_from_coordinates(x, y, t)
    }

    fn next_neighbour_index(&self, index: usize, direction: Direction) -> usize {
        let (mut x, mut y, mut t) = Self::get_coordinates_from_index(index);
        match direction {
            Direction::X => x = (x + 1) % MAX_X,
            Direction::Y => y = (y + 1) % MAX_Y,
            Direction::T => t = (t + 1) % MAX_T,
        }

        Self::get_index_from_coordinates(x, y, t)
    }

    fn get_index_from_coordinates(x: usize, y: usize, t: usize) -> usize {
        MAX_X * MAX_Y * t + MAX_X * y + x
    }

    fn get_coordinates_from_index(index: usize) -> (usize, usize, usize) {
        (
            Self::get_x_from_index(index),
            Self::get_y_from_index(index),
            Self::get_t_from_index(index),
        )
    }

    fn get_x_from_index(index: usize) -> usize {
        index % MAX_X
    }

    fn get_y_from_index(index: usize) -> usize {
        (index / MAX_X) % MAX_Y
    }

    fn get_t_from_index(index: usize) -> usize {
        index / (MAX_X * MAX_Y) % MAX_T
    }
}

pub trait LatticeValues<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    T: std::fmt::Debug,
    Self: Lattice<MAX_X, MAX_Y, MAX_T>,
{
    fn prev_neighbour_value(&self, index: usize, direction: Direction) -> &T {
        let neighbour_index = Self::prev_neighbour_index(&self, index, direction);
        Self::get_value(&self, neighbour_index)
    }

    fn next_neighbour_value(&self, index: usize, direction: Direction) -> &T {
        let neighbour_index = Self::next_neighbour_index(&self, index, direction);
        Self::get_value(&self, neighbour_index)
    }

    fn get_value(&self, index: usize) -> &T;

    fn print_values(&self);

    fn print_values_formated(&self) {
        for t in 0..MAX_T {
            println!("t = {}", t);
            for y in 0..MAX_Y {
                print!("[");
                for x in 0..MAX_X {
                    if x == MAX_X - 1 {
                        println!(
                            "{:?} ]",
                            self.get_value(Self::get_index_from_coordinates(x, y, t))
                        );
                    } else {
                        print!(
                            "{:?}, ",
                            self.get_value(Self::get_index_from_coordinates(x, y, t))
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn test_index_coordinates_conversion() {
    let (x, y, t) = Lattice3d::<i32, 4, 5, 3>::get_coordinates_from_index(29);
    let index = Lattice3d::<i32, 4, 5, 3>::get_index_from_coordinates(x, y, t);

    assert_eq!(x, 1);
    assert_eq!(y, 2);
    assert_eq!(t, 1);
    assert_eq!(index, 29);
}

#[test]
fn test_neighbour_index() {
    let lattice = Lattice3d::<i32, 4, 5, 3>::default();
    let center: usize = 19;

    assert_eq!(lattice.next_neighbour_index(center, Direction::X), 16);
    assert_eq!(lattice.next_neighbour_index(center, Direction::Y), 3);
    assert_eq!(lattice.next_neighbour_index(center, Direction::T), 39);
    assert_eq!(lattice.prev_neighbour_index(center, Direction::X), 18);
    assert_eq!(lattice.prev_neighbour_index(center, Direction::Y), 15);
    assert_eq!(lattice.prev_neighbour_index(center, Direction::T), 59);
}

#[test]
fn test_lattice3d_default() {
    let lattice = Lattice3d::<i32, 4, 5, 3>::default();
    assert_eq!(lattice.0, [0_i32; 60]);
}
