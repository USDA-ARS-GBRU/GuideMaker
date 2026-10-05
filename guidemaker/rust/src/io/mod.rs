pub mod fasta;
pub mod gff;

use anyhow::{Context, Result};
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use flate2::read::GzDecoder;
use zstd::stream::Decoder as ZstdDecoder;

/// Safely extract column values as a vector of String across Categorical, String, or Enum types
pub fn column_to_string_vec(df: &DataFrame, name: &str) -> Result<Vec<String>> {
    let col = df.column(name)?;
    let casted = col.cast(&DataType::String)?;
    let str_ca = casted.str()?;
    Ok(str_ca
        .into_iter()
        .map(|opt| opt.unwrap_or("").to_string())
        .collect())
}

/// Transparently open a file with optional gzip or zstd decompression
pub fn open_compressed_reader(path: &Path) -> Result<Box<dyn Read + Send>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open sequence file at {:?}", path))?;
    let mut buf_reader = BufReader::new(file);

    let mut magic = [0u8; 4];
    let n = buf_reader.read(&mut magic)?;

    let combined_reader = std::io::Cursor::new(magic[..n].to_vec()).chain(buf_reader);

    if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        let gz_decoder = GzDecoder::new(combined_reader);
        Ok(Box::new(gz_decoder))
    } else if n >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd {
        let zstd_decoder = ZstdDecoder::new(combined_reader)
            .context("Failed to initialize Zstd decoder")?;
        Ok(Box::new(zstd_decoder))
    } else {
        Ok(Box::new(combined_reader))
    }
}

/// Write Polars DataFrame to CSV
pub fn write_csv(df: &mut DataFrame, path: &Path) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create CSV file at {:?}", path))?;
    CsvWriter::new(file).finish(df)?;
    Ok(())
}

/// Write Polars DataFrame to Parquet
pub fn write_parquet(df: &mut DataFrame, path: &Path) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create Parquet file at {:?}", path))?;
    ParquetWriter::new(file).finish(df)?;
    Ok(())
}
