//! This module implements a data structure to save the field configuration
//! data of the Markov chain Monte Carlo in order to calculate observables
//! and create plots.
//!
//! It implements a

use crate::{
    fields::{Action, BondsFieldNew, Field, HeightVariable},
    kahan::KahanSummation,
    lattice::Lattice,
};

pub trait ObservableField<const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    fn new(lattice: &Lattice<D, SIZE>) -> Self;
    fn update<T: HeightVariable<T>>(&mut self, field: &Field<T, D, SIZE>);
}

#[derive(Debug)]
pub struct EnergyObservableField<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    bonds: BondsFieldNew<'a, KahanSummation<f64>, D, SIZE>,
}

impl<'a, const D: usize, const SIZE: usize> ObservableField<D, SIZE>
    for EnergyObservableField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn new(lattice: &Lattice<D, SIZE>) -> Self {
        EnergyObservableField {
            bonds: BondsFieldNew::new(lattice),
        }
    }

    fn update<T: HeightVariable<T>>(&mut self, field: &Field<T, D, SIZE>) {
        for index in 0..SIZE {
            for direction in 0..D {
                let action = field.bond_action(index, direction);
                self.bonds
                    .get_value_mut(index, direction)
                    .add(action.into());
            }
        }
    }
}

#[derive(Debug)]
pub struct DifferenceObservableField<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    bonds: BondsFieldNew<'a, KahanSummation<f64>, D, SIZE>,
}

impl<'a, const D: usize, const SIZE: usize> ObservableField<D, SIZE>
    for DifferenceObservableField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn new(lattice: &Lattice<D, SIZE>) -> Self {
        DifferenceObservableField {
            bonds: BondsFieldNew::new(lattice),
        }
    }

    fn update<T: HeightVariable<T>>(&mut self, field: &Field<T, D, SIZE>) {
        for index in 0..SIZE {
            for direction in 0..D {
                let diff = field.bond_difference(index, direction);
                self.bonds.get_value_mut(index, direction).add(diff.into());
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
