use serde::{Deserialize, Serialize};

use crate::{
    field::{Field, HeightVariable},
    wilson::Wilsonloop,
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum ActionType {
    StandardAction,
    WilsonAction(Wilsonloop),
}

impl<T, const D: usize, const SIZE: usize> Action<T, D, SIZE> for ActionType
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn bond_formula(&self, site1_val: T, site2_val: T) -> T {
        match self {
            ActionType::StandardAction => <StandardAction as Action<T, D, SIZE>>::bond_formula(
                &StandardAction,
                site1_val,
                site2_val,
            ),
            ActionType::WilsonAction(wilson) => <WilsonAction as Action<T, D, SIZE>>::bond_formula(
                &WilsonAction(wilson),
                site1_val,
                site2_val,
            ),
        }
    }
}

pub struct StandardAction;
pub struct WilsonAction<'a>(&'a Wilsonloop);

pub trait Action<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn bond_formula(&self, site1_val: T, site2_val: T) -> T;

    /// The standard action is determined as the difference squared.
    fn bond_action(&self, field: &Field<T, D, SIZE>, site1: usize, site2: usize) -> T {
        Self::bond_formula(self, field.values[site1], field.values[site2])
    }

    fn direction_bond_formula(
        &self,
        field: &Field<T, D, SIZE>,
        value: T,
        index: usize,
        direction: usize,
    ) -> T {
        let neighbour: T = field.values[field.lattice.get_neighbours_array(index)[direction]];
        self.bond_formula(value, neighbour)
    }

    fn direction_bond_action(
        &self,
        field: &Field<T, D, SIZE>,
        index: usize,
        direction: usize,
    ) -> T {
        Self::direction_bond_formula(&self, field, field.values[index], index, direction)
    }

    fn calculate_assumed_action(&self, field: &Field<T, D, SIZE>, index: usize, value: T) -> T {
        let mut sum: T = T::default();
        for neighbour in field.lattice.get_neighbours_array(index) {
            sum = sum + Self::bond_formula(self, value, field.values[neighbour]);
        }
        sum
    }

    fn integer_action(&self, field: &Field<T, D, SIZE>) -> T {
        let mut sum: T = T::default();
        for index in 0..SIZE {
            for neighbour in field.lattice.pos_neighbours_array(index) {
                sum = sum + self.bond_action(field, index, neighbour);
            }
        }
        sum
    }

    fn action_observable(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64 {
        Into::<f64>::into(self.integer_action(field)) * temp
    }
}

impl<T, const D: usize, const SIZE: usize> Action<T, D, SIZE> for StandardAction
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn bond_formula(&self, site1_val: T, site2_val: T) -> T {
        let diff: T = site1_val - site2_val;
        diff * diff
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Action<T, D, SIZE> for WilsonAction<'a>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    /// The wilson action is the standard action with a offset if the bond
    /// passes the area enclosed by the wilson loop.
    fn bond_formula(&self, _: T, _: T) -> T {
        panic!("No bond formula for wilson action only dependant on the values");
    }

    fn bond_action(&self, _: &Field<T, D, SIZE>, _: usize, _: usize) -> T {
        panic!("No bond wilson action without direction");
    }

    fn direction_bond_formula(
        &self,
        field: &Field<T, D, SIZE>,
        value: T,
        index: usize,
        direction: usize,
    ) -> T {
        self.0
            .direction_bond_formula(field, value, index, direction)
    }

    fn direction_bond_action(
        &self,
        field: &Field<T, D, SIZE>,
        index: usize,
        direction: usize,
    ) -> T {
        self.0.direction_bond_action(field, index, direction)
    }

    fn calculate_assumed_action(&self, field: &Field<T, D, SIZE>, index: usize, value: T) -> T {
        self.0.calculate_assumed_action(field, index, value)
    }

    fn integer_action(&self, field: &Field<T, D, SIZE>) -> T {
        self.0.integer_action(field)
    }
}
