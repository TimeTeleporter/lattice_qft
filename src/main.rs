#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use forms3d::ZeroForm;
use lattice3d::Lattice3d;

const MAX_X: usize = 10;
const MAX_Y: usize = 10;
const MAX_T: usize = 10;

mod forms3d;
mod lattice3d;

fn main() {
    let lattice = Lattice3d::<MAX_X, MAX_Y, MAX_T>::default();

    let mut m_field = ZeroForm::<i32, MAX_X, MAX_Y, MAX_T>::from(&lattice);

    m_field.set_value(1, 3, 1, 2);

    m_field.set_value(2, 5, 7, 1);

    m_field.print_values_formated();

    m_field.set_value(3, 5, 7, 1);

    let k_field = ZeroForm::<i32, MAX_X, MAX_Y, MAX_T>::from(&lattice);
}
