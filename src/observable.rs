use serde::{Deserialize, Serialize};

use crate::{
    action::{Action, ActionType},
    field::{Field, HeightVariable},
    wilson::Wilsonloop,
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum ObservableType {
    ActionObservable(ActionType),
    WilsonObservable(Wilsonloop),
    SizeNormalized(ActionType),
    BondNormalized(ActionType),
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

pub trait Observable<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn observe(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64;
}

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
