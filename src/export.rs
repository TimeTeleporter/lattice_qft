use std::{error::Error, fs::File};

use serde::{Deserialize, Serialize};

use crate::computation::{ComputationSummary, FieldExport3d};

pub trait CsvData
where
    Self: Serialize + for<'a> Deserialize<'a>,
{
    fn read_write_csv(self, path: &str, has_headers: bool) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<Self> = Self::fetch_csv_data(path, has_headers)?;

        // Append the new entry
        storage.push(self);

        // Write to the file
        write_to_csv(path, storage)?;

        Ok(())
    }

    fn overwrite_csv(self, path: &str) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<Self> = Vec::new();

        // Append the new entry
        storage.push(self);

        write_to_csv(path, storage)?;

        Ok(())
    }

    fn fetch_csv_data(path: &str, has_headers: bool) -> Result<Vec<Self>, Box<dyn Error>> {
        let mut storage: Vec<Self> = Vec::new();

        read_from_csv(path, &mut storage, has_headers)?;

        Ok(storage)
    }
}

fn write_to_csv<T>(path: &str, storage: Vec<T>) -> Result<(), Box<dyn Error>>
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    let mut wtr = csv::Writer::from_path(path)?;
    for data in storage {
        wtr.serialize(data)?;
    }
    wtr.flush()?;
    Ok(())
}

fn read_from_csv<T>(
    path: &str,
    storage: &mut Vec<T>,
    has_headers: bool,
) -> Result<(), Box<dyn Error>>
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    let file: File = match File::open(path) {
        Ok(file) => file,
        Err(err) => match err.kind() {
            std::io::ErrorKind::NotFound => {
                File::create(path)?;
                File::open(path)?
            }
            _ => Err(err)?,
        },
    };
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(has_headers)
        .from_reader(file);
    Ok(for result in rdr.deserialize() {
        let test_data: T = result?;
        storage.push(test_data);
    })
}

pub fn clean_csv(path: &str) -> Result<(), Box<dyn Error>> {
    let storage: Vec<()> = Vec::new();
    write_to_csv(path, storage)?;
    Ok(())
}

pub fn get_correlation_fn(index: usize) -> Result<Vec<f64>, Box<dyn Error>> {
    let path: &str =
        &(crate::PLOT_PATH_INCOMPLETE.to_owned() + &"correlation_" + &index.to_string() + &".csv");
    f64::fetch_csv_data(path, false)
}

impl CsvData for ComputationSummary {}
impl CsvData for FieldExport3d<f64> {}
impl CsvData for f64 {}

impl<T: CsvData> CsvData for Vec<T> {
    fn read_write_csv(mut self, path: &str, has_headers: bool) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<T> = <T as CsvData>::fetch_csv_data(path, has_headers)?;

        // Append the new entry
        storage.append(&mut self);

        // Write to the file
        write_to_csv(path, storage)?;

        Ok(())
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FitResult {
    index: usize,
    m: f64,
    n: f64,
    a: f64,
    b: f64,
}

impl FitResult {
    pub fn new(index: usize, m: f64, n: f64, a: f64, b: f64) -> FitResult {
        FitResult { index, m, n, a, b }
    }
}

impl CsvData for FitResult {}

#[derive(Debug, Serialize, Deserialize)]
pub struct CorrelationLengths {
    index: usize,
    m12: f64,
    m23: f64,
    m13: f64,
    corr12: f64,
    corr23: f64,
    corr13: f64,
}

impl CorrelationLengths {
    pub fn new(index: usize, m12: f64, m23: f64, m13: f64) -> CorrelationLengths {
        CorrelationLengths {
            index,
            m12,
            m23,
            m13,
            corr12: 1.0 / m12,
            corr23: 1.0 / m23,
            corr13: 1.0 / m13,
        }
    }
}

impl CsvData for CorrelationLengths {}
