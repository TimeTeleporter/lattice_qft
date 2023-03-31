use std::{error::Error, fs::File};

use serde::{Deserialize, Serialize};

use crate::computation::{ComputationSummary, FieldExport3d};

pub trait CsvData
where
    Self: Serialize + for<'a> Deserialize<'a>,
{
    fn read_write_csv(self, path: &str) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<Self> = Self::fetch_csv_data(path)?;

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

    fn fetch_csv_data(path: &str) -> Result<Vec<Self>, Box<dyn Error>> {
        let mut storage: Vec<Self> = Vec::new();

        read_from_csv(path, &mut storage)?;

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

fn read_from_csv<T>(path: &str, storage: &mut Vec<T>) -> Result<(), Box<dyn Error>>
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
    let mut rdr = csv::Reader::from_reader(file);
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

impl CsvData for ComputationSummary {}
impl CsvData for FieldExport3d<f64> {}
impl CsvData for Vec<FieldExport3d<f64>> {
    fn read_write_csv(mut self, path: &str) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<FieldExport3d<f64>> =
            <FieldExport3d<f64> as CsvData>::fetch_csv_data(path)?;

        // Append the new entry
        storage.append(&mut self);

        // Write to the file
        write_to_csv(path, storage)?;

        Ok(())
    }
}
