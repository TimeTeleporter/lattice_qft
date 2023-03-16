//! This module implements a data structure to save the field configuration
//! data of the Markov chain Monte Carlo in order to calculate observables
//! and create plots.

use crate::{
    fields::{BondsField, Field},
    kahan::KahanSummation,
    lattice::Lattice,
};

pub trait ObservableField<const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    fn new(lattice: &Lattice<D, SIZE>) -> Self;
    fn update<T: Into<f64>>(&mut self, field: &Field<T, D, SIZE>);
}

#[derive(Debug)]
pub struct EnergyObservableField<'a, const D: usize, const SIZE: usize>(
    BondsField<'a, KahanSummation<f64>, D, SIZE>,
)
where
    [(); D * 2_usize]:;

impl<'a, const D: usize, const SIZE: usize> ObservableField<D, SIZE>
    for EnergyObservableField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn new(lattice: &Lattice<D, SIZE>) -> Self {
        EnergyObservableField(BondsField::new(lattice))
    }

    fn update<T: Into<f64>>(&mut self, field: &Field<T, D, SIZE>) {
        for index in 0..SIZE {
            for direction in 0..D {
                self.0.set_value(index, direction, value)
            }
        }
    }
}

#[test]
fn test_observable_bond_field() {
    const D: usize = 3;
    const SIZE: usize = 4 * 5 * 6;
    let lattice: Lattice<D, SIZE> = Lattice::new([4, 5, 6]);
    let diff_field: EnergyObservableField<D, SIZE> = EnergyObservableField::new(&lattice);
    let energy_field: EnergyObservableField<D, SIZE> = EnergyObservableField::new(&lattice);
}
