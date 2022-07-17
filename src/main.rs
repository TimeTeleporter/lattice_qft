#![feature(generic_const_exprs)]

// Number of spacetime dimensions of our lattice
pub const DIMENSIONS: usize = 3;

const CUBE_SIDE: usize = 5;

// The size of the spacetime volume
pub const MAX_X: usize = CUBE_SIDE;
pub const MAX_Y: usize = CUBE_SIDE;
pub const MAX_T: usize = CUBE_SIDE;

mod form;
mod lattice3d;

use form::Form;
use lattice3d::{Direction, Lattice3d};

type MyInt = i32;

type ZeroForm = Form<MyInt, DIMENSIONS, 0>;
type OneForm = Form<MyInt, DIMENSIONS, 1>;
type TwoForm = Form<MyInt, DIMENSIONS, 2>;
type ThreeForm = Form<MyInt, DIMENSIONS, 3>;

fn main() {
    let mut m_field = Lattice3d::<OneForm, MAX_X, MAX_Y, MAX_T>::new();

    for (x, x_slice) in m_field.0.iter_mut().enumerate() {
        for (y, xy_slice) in x_slice.iter_mut().enumerate() {
            for (t, mut value) in xy_slice.iter_mut().enumerate() {
                let mut count: usize = 0;
                let ary = [x, y, t];
                while count < 3 {
                    value.components[count] = ary[count] as MyInt;
                    count += 1;
                }
            }
        }
    }

    println!("{:?}", m_field.value(2, 2, 2));
    println!("{:?}", m_field.prev_neighbour_value(2, 2, 2, Direction::X));
    println!("{:?}", m_field.next_neighbour_value(2, 2, 2, Direction::T));
    println!("{:?}", m_field.next_neighbour_value(2, 2, 2, Direction::Y));
}
