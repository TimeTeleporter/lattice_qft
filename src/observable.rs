use std::fmt::Display;

use serde::{Deserialize, Serialize};

use crate::{
    field::Field,
    heightfield::{Action, HeightVariable},
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum ObservableType {
    Action,
    SizeNormalizedAction,
    BondNormalizedAction,
}

impl Display for ObservableType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObservableType::Action => write!(f, "action                "),
            ObservableType::SizeNormalizedAction => write!(f, "action size normalized"),
            ObservableType::BondNormalizedAction => write!(f, "action bond normalized"),
        }
    }
}

/// This trait implements methods to get a observable from a [HeightField].
pub trait Observable<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn observe(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64;
}

impl<T, const D: usize, const SIZE: usize> Observable<T, D, SIZE> for ObservableType
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn observe(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64 {
        match self {
            ObservableType::Action => field.action_observable(temp),
            ObservableType::SizeNormalizedAction => {
                field.action_observable(temp) / f64::from(SIZE as u32)
            }
            ObservableType::BondNormalizedAction => {
                field.action_observable(temp) / f64::from(SIZE as u32 * 3)
            }
        }
    }
}
