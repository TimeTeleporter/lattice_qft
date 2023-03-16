//! This module implements Kahan summation according to [https://en.wikipedia.org/wiki/Kahan_summation_algorithm].
use std::ops::{Add, Sub};

/// A data structure using Kahan summation of floating point numbers to update a running average.
#[derive(Debug, Clone, Default)]
pub struct KahanSummation<T> {
    value: T,
    kahan: T,
    entries: usize,
}

impl<T> KahanSummation<T> {
    pub fn new() -> Self
    where
        T: Default,
    {
        KahanSummation::default()
    }

    /// Adds a value according to Kahan summation
    pub fn add(&mut self, value: T)
    where
        T: Add<Output = T> + Sub<Output = T> + Clone,
    {
        let y: T = value - self.kahan.clone();
        let t: T = self.value.clone() + y.clone();
        self.kahan = (t.clone() - self.value.clone()) - y;
        self.value = t;
        self.entries = self.entries + 1;
    }

    /// Returns the result of the Kahan summation
    pub fn sum(&self) -> T
    where
        T: Clone,
    {
        self.value.clone()
    }

    /// Returns the average of the added values
    pub fn mean(&self) -> f64
    where
        T: Into<f64> + Clone,
    {
        self.value.clone().into() / self.entries as f64
    }
}

#[test]
fn test_kahan_default() {
    let obs: KahanSummation<f64> = KahanSummation::default();
    assert_eq!(obs.value, 0.0);
    assert_eq!(obs.kahan, 0.0);
    assert_eq!(obs.entries, 0);
}

#[test]
fn test_kahan() {
    let mut k: KahanSummation<f32> = KahanSummation::new();
    k.add(3.1);
    k.add(-3.1);
    assert_eq!(0.0, k.sum());
}
