use std::{
    fmt::Display,
    ops::{Add, Sub},
};

use crate::{
    fields::{Action, BondsField, Field, HeightVariable},
    lattice::Lattice,
};

// - ObservableValue ----------------------------------------------------------

/// A data structure using kahan summation to update a running average.
#[derive(Debug, Clone, Default)]
pub struct ObservableValue<T> {
    value: T,
    kahan: T,
    entries: usize,
}

impl<T> ObservableValue<T> {
    pub fn update(&mut self, value: T)
    where
        T: Add<Output = T> + Sub<Output = T> + Clone,
    {
        let y: T = value - self.kahan.clone();
        let t: T = self.value.clone() + y.clone();
        self.kahan = (t.clone() - self.value.clone()) - y;
        self.value = t;
        self.entries = self.entries + 1;
    }

    pub fn result(&self) -> f64
    where
        T: Into<f64> + Clone,
    {
        self.value.clone().into() / self.entries as f64
    }
}

// - ObservableType -----------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub enum ObservableType {
    Action,
    SizeNormalizedAction,
    BondNormalizedAction,
}

impl ObservableType {
    fn bond_observe<T: HeightVariable<T>, F: Action<T>>(
        &self,
        field: &F,
        index: usize,
        direction: usize,
    ) -> T {
        match self {
            ObservableType::Action
            | ObservableType::SizeNormalizedAction
            | ObservableType::BondNormalizedAction => field.bond_action(index, direction),
        }
    }
}

impl ObservableType {
    pub fn observe<'a, T, const D: usize, const SIZE: usize>(
        &self,
        field: &'a Field<T, D, SIZE>,
        temp: f64,
    ) -> f64
    where
        T: HeightVariable<T>,
        [(); D * 2_usize]:,
    {
        field.action_observable(temp)
            * match self {
                ObservableType::Action => 1.0,
                ObservableType::SizeNormalizedAction => 1.0 / SIZE as f64,
                ObservableType::BondNormalizedAction => 1.0 / (SIZE * D) as f64,
            }
    }
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

// - Observable ---------------------------------------------------------------

pub struct Observable<'a, T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    values: BondsField<'a, ObservableValue<T>, D, SIZE>,
    obstype: ObservableType,
    temp: f64,
}

impl<'a, T, const D: usize, const SIZE: usize> Observable<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    pub fn new(
        lattice: &'a Lattice<D, SIZE>,
        obstype: ObservableType,
        temp: f64,
    ) -> Observable<'a, T, D, SIZE> {
        Observable {
            values: BondsField::<'a, ObservableValue<T>, D, SIZE>::new(lattice),
            obstype,
            temp,
        }
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Observable<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    pub fn update(&mut self, field: &Field<T, D, SIZE>) {
        for index in 0..SIZE {
            for direction in 0..D {
                let value: T = self.obstype.bond_observe(field, index, direction);
                self.values.values[index][direction].update(value);
            }
        }
    }

    pub fn result(&self) -> f64 {
        let mut res: Vec<f64> = Vec::with_capacity(SIZE * D);
        for index in 0..SIZE {
            for direction in 0..D {
                res.push(self.values.values[index][direction].result());
            }
        }
        let res: f64 = res.into_iter().sum();
        let res: f64 = res * self.temp;
        match self.obstype {
            ObservableType::Action => res,
            ObservableType::SizeNormalizedAction => res / (SIZE as f64),
            ObservableType::BondNormalizedAction => res / ((D * SIZE) as f64),
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

#[test]
fn requirements_observable() {
    use crate::lattice::Lattice;

    const D: usize = 3;
    const SIZE_ARY: [usize; 3] = [3, 3, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];
    const TEMP: f64 = 0.1;

    let lattice: Lattice<D, SIZE> = Lattice::new([3, 3, 3]);
    let field: Field<i32, D, SIZE> = Field::new(&lattice);

    let obs_type: ObservableType = ObservableType::Action;
    let mut observable: Observable<i32, D, SIZE> = Observable::new(&lattice, obs_type, TEMP);

    observable.update(&field);

    let _: f64 = observable.result();
}
