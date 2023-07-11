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
    pub fn get_sum(&self) -> T
    where
        T: Clone,
    {
        self.value.clone()
    }

    pub fn get_mean(&self) -> f64
    where
        T: Clone + Into<f64>,
    {
        match self.entries {
            0 => 0.0,
            entries => self.value.clone().into() / (entries as f64),
        }
    }
}

impl Into<f64> for KahanSummation<f64> {
    fn into(self) -> f64 {
        self.get_mean()
    }
}

#[test]
fn test_kahan_new() {
    let mut obs: KahanSummation<f64> = KahanSummation::new();
    assert_eq!(obs.get_sum(), 0.0);
    assert_eq!(obs.get_mean(), 0.0);
    assert_eq!(obs.kahan, 0.0);
    obs.add(1.0);
    assert_eq!(obs.get_sum(), 1.0);
    assert_eq!(obs.get_mean(), 1.0);
    assert_eq!(obs.kahan, 0.0);
}

#[test]
fn test_kahan() {
    let mut k: KahanSummation<f32> = KahanSummation::new();
    k.add(3.1);
    k.add(-3.1);
    assert_eq!(0.0, k.get_sum());
}

#[derive(Debug, Default, Clone)]
pub struct WelfordsAlgorithm64 {
    sum: KahanSummation<f64>,
    diff_square_sum: KahanSummation<f64>,
    entries: usize,
}

impl WelfordsAlgorithm64 {
    pub fn new() -> Self {
        WelfordsAlgorithm64::default()
    }

    pub fn get_mean(&self) -> f64 {
        match self.entries {
            0 => 0.0,
            entries => self.sum.get_sum() / (entries as f64),
        }
    }

    pub fn get_count(&self) -> usize {
        self.entries
    }

    pub fn update(&mut self, value: f64) {
        let old_mean: f64 = self.get_mean();
        self.sum.add(value);
        self.entries = self.entries + 1;
        let new_mean: f64 = self.get_mean();
        self.diff_square_sum
            .add((value - old_mean) * (value - new_mean));
    }

    /// This method yields the biased sample variance. For the unbiased sample
    /// variance use [get_var_unbiased].
    pub fn get_var(&self) -> f64 {
        match self.entries {
            0 => 0.0,
            entries => self.diff_square_sum.get_sum() / (entries as f64),
        }
    }

    pub fn get_sd(&self) -> f64 {
        self.get_var().sqrt()
    }

    /// This method yields the unbiased sample variance.
    pub fn get_var_unbiased(&self) -> f64 {
        match self.entries {
            0 | 1 => 0.0,
            entries => self.diff_square_sum.get_sum() / ((entries - 1) as f64),
        }
    }

    pub fn get_sd_unbiased(&self) -> f64 {
        self.get_var_unbiased().sqrt()
    }

    /// This method yield the standard error of the mean.
    pub fn get_standard_error(&self) -> f64 {
        self.get_var_unbiased() / (self.entries as f64).sqrt()
    }
}

#[test]
fn test_welfords_mean() {
    let values: Vec<f64> = vec![1.0, 0.0, 1.0, 2.0];
    let mean: f64 = values.iter().cloned().sum::<f64>() / values.len() as f64;
    let mut welford: WelfordsAlgorithm64 = WelfordsAlgorithm64::new();
    values
        .iter()
        .cloned()
        .for_each(|value| welford.update(value));
    assert_eq!(mean, welford.get_mean());
}

#[test]
fn test_welfords_sd() {
    let values: Vec<f64> = vec![1.0, 0.0, 1.0, 2.0, 25.0];
    let mean: f64 = values.iter().cloned().sum::<f64>() / values.len() as f64;
    let var: f64 = values
        .iter()
        .cloned()
        .map(|x| (x - mean) * (x - mean))
        .sum::<f64>()
        / values.len() as f64;
    let mut welford: WelfordsAlgorithm64 = WelfordsAlgorithm64::new();
    values.iter().cloned().for_each(|value| {
        welford.update(value);
    });
    assert_eq!(var, welford.get_var());
    assert_eq!(var.sqrt(), welford.get_sd());
}
