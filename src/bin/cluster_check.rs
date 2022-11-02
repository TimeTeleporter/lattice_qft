#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use lattice_qft::{cluster::Cluster, field::Field, lattice::Lattice, pause};

const D: usize = 2;
const SIZE_ARY: [usize; D] = [3, 3];
const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1];

const TEMP: f64 = 0.1;

fn main() {
    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let field: Field<i8, D, SIZE> = Field::random(&lattice);
    let mut field: Field<i32, D, SIZE> = Field::from_field(field);

    field.print_values_formated(SIZE_ARY);
    pause();
    field.cluster_sweep(TEMP);
    field.print_values_formated(SIZE_ARY);
}
