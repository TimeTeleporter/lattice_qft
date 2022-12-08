use std::fmt::Display;

use serde::{Deserialize, Serialize};

use crate::heightfield::HeightVariable;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum ObservableType {
    ActionObservable(ActionType),
    WilsonObservable(Wilsonloop),
    SizeNormalized(ActionType),
    BondNormalized(ActionType),
}

impl Display for ObservableType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObservableType::ActionObservable(action) => write!(f, "{} observable", action),
            ObservableType::WilsonObservable(wilson) => write!(f, "{} observable", wilson),
            ObservableType::SizeNormalized(action) => {
                write!(f, "{} size normalized observable", action)
            }
            ObservableType::BondNormalized(action) => {
                write!(f, "{} bond normalized observable", action)
            }
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
            ObservableType::ActionObservable(action) => {
                ActionObservable(action).observe(field, temp)
            }
            ObservableType::WilsonObservable(wilson) => {
                WilsonObservable(wilson).observe(field, temp)
            }
            ObservableType::SizeNormalized(action) => {
                ActionObservable(action).observe(field, temp) / SIZE as f64
            }
            ObservableType::BondNormalized(action) => {
                ActionObservable(action).observe(field, temp) / (SIZE * 3) as f64
            }
        }
    }
}

struct ActionObservable<'a>(&'a ActionType);
struct WilsonObservable<'a>(&'a Wilsonloop);

impl<'a, T, const D: usize, const SIZE: usize> Observable<T, D, SIZE> for ActionObservable<'a>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn observe(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64 {
        self.0.action_observable(field, temp)
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Observable<T, D, SIZE> for WilsonObservable<'a>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn observe(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64 {
        self.0.action_observable(field, temp)
    }
}
