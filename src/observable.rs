use std::{
    fmt::Display,
    ops::{Add, Sub},
};

use crate::{
    field::Field,
    heightfield::{Action, HeightVariable},
};

#[derive(Debug, Clone, Default)]
pub struct ObservableValue<T: Default + Clone> {
    value: T,
    kahan: T,
    entries: usize,
}

impl<T: Default + Add<Output = T> + Sub<Output = T> + Clone> ObservableValue<T> {
    pub fn update(&mut self, value: T) {
        let y: T = value - self.kahan.clone();
        let t: T = self.value.clone() + y.clone();
        self.kahan = (t.clone() - self.value.clone()) - y;
        self.value = t;
        self.entries = self.entries + 1;
    }

    pub fn result(&self) -> f64
    where
        T: Into<f64>,
    {
        self.value.clone().into() / self.entries as f64
    }
}

#[derive(Debug, Clone)]
pub enum ObservableType {
    Action(ObservableValue<f64>),
    SizeNormalizedAction(ObservableValue<f64>),
    BondNormalizedAction(ObservableValue<f64>),
}

impl Display for ObservableType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObservableType::Action(_) => write!(f, "action                "),
            ObservableType::SizeNormalizedAction(_) => write!(f, "action size normalized"),
            ObservableType::BondNormalizedAction(_) => write!(f, "action bond normalized"),
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
    fn update(&mut self, field: &Field<T, D, SIZE>, temp: f64);
    fn result(&self) -> f64;
}

impl<T, const D: usize, const SIZE: usize> Observable<T, D, SIZE> for ObservableType
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn observe(&self, field: &Field<T, D, SIZE>, temp: f64) -> f64 {
        match self {
            ObservableType::Action(_) => field.action_observable(temp),
            ObservableType::SizeNormalizedAction(_) => {
                field.action_observable(temp) / f64::from(SIZE as u32)
            }
            ObservableType::BondNormalizedAction(_) => {
                field.action_observable(temp) / f64::from(SIZE as u32 * 3)
            }
        }
    }

    fn update(&mut self, field: &Field<T, D, SIZE>, temp: f64) {
        let value: f64 = self.observe(field, temp);
        match self {
            ObservableType::Action(obs)
            | ObservableType::SizeNormalizedAction(obs)
            | ObservableType::BondNormalizedAction(obs) => obs.update(value),
        }
    }

    fn result(&self) -> f64 {
        match self {
            ObservableType::Action(value)
            | ObservableType::SizeNormalizedAction(value)
            | ObservableType::BondNormalizedAction(value) => value.result(),
        }
    }
}

#[test]
fn test_observable_value_default() {
    let obs: ObservableValue<f64> = ObservableValue::default();
    assert_eq!(obs.value, 0.0);
    assert_eq!(obs.kahan, 0.0);
    assert_eq!(obs.entries, 0);
}
