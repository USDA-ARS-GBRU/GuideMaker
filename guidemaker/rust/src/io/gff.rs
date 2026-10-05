use anyhow::{Context, Result};
use polars::prelude::*;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Representation of a genomic feature record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureRecord {
    pub chrom: String,
    pub feature_start: u32,
    pub feature_end: u32,
    pub strand: bool,
    pub feature_id: String,
    pub feature_type: String,
}

/// Parse GFF/GTF features from reader
pub fn parse_gff_gtf_features<R: Read>(reader: R) -> Result<Vec<FeatureRecord>> {
    let buf_reader = BufReader::new(reader);
    let mut features = Vec::new();

    for line_res in buf_reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let chrom = cols[0].to_string();
        let feature_type = cols[2].to_string();
        let start_1based: u32 = match cols[3].parse() {
            Ok(val) => val,
            Err(_) => continue,
        };
        let end_1based: u32 = match cols[4].parse() {
            Ok(val) => val,
            Err(_) => continue,
        };

        let feature_start = start_1based.saturating_sub(1);
        let feature_end = end_1based;
        let strand = cols[6] != "-";

        let attrs = cols[8];
        let feature_id = extract_gff_gtf_id(attrs, &chrom, feature_start, feature_end, &feature_type);

        features.push(FeatureRecord {
            chrom,
            feature_start,
            feature_end,
            strand,
            feature_id,
            feature_type,
        });
    }

    Ok(features)
}

fn extract_gff_gtf_id(
    attrs: &str,
    chrom: &str,
    start: u32,
    end: u32,
    ftype: &str,
) -> String {
    let mut found_id = String::new();

    for part in attrs.split(';') {
        let trimmed = part.trim();
        if let Some((key, val)) = trimmed.split_once('=') {
            let k = key.trim();
            let clean_val = val.trim_matches('"').trim();
            if !clean_val.is_empty() {
                if k.eq_ignore_ascii_case("locus_tag") || k.eq_ignore_ascii_case("ID") {
                    return clean_val.to_string();
                } else if (k.eq_ignore_ascii_case("gene_id") || k.eq_ignore_ascii_case("Name") || k.eq_ignore_ascii_case("transcript_id"))
                    && found_id.is_empty()
                {
                    found_id = clean_val.to_string();
                }
            }
        }
    }

    if !found_id.is_empty() {
        return found_id;
    }

    for part in attrs.split(';') {
        let trimmed = part.trim();
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() >= 2 {
            let k = tokens[0];
            let clean_val = tokens[1].trim_matches('"').trim();
            if !clean_val.is_empty() {
                if k.eq_ignore_ascii_case("locus_tag") || k.eq_ignore_ascii_case("gene_id") {
                    return clean_val.to_string();
                } else if k.eq_ignore_ascii_case("transcript_id") && found_id.is_empty() {
                    found_id = clean_val.to_string();
                }
            }
        }
    }

    if !found_id.is_empty() {
        return found_id;
    }

    format!("{}_{}_{}_{}", chrom, start, end, ftype)
}

/// Parse GenBank FEATURES section into FeatureRecord list
pub fn parse_genbank_features<R: Read>(reader: R) -> Result<Vec<FeatureRecord>> {
    let buf_reader = BufReader::new(reader);
    let mut features = Vec::new();

    let mut current_chrom = "unknown_contig".to_string();
    let mut in_features = false;
    let mut current_type = String::new();
    let mut current_loc = String::new();
    let mut current_id = String::new();

    for line_res in buf_reader.lines() {
        let line = line_res?;

        if line.starts_with("LOCUS") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() > 1 {
                current_chrom = parts[1].to_string();
            }
            in_features = false;
        } else if line.starts_with("FEATURES") {
            in_features = true;
        } else if line.starts_with("ORIGIN") || line.starts_with("//") {
            if !current_type.is_empty() {
                if let Some((f_start, f_end, strand)) = parse_genbank_location(&current_loc) {
                    let fid = if !current_id.is_empty() {
                        std::mem::take(&mut current_id)
                    } else {
                        format!("{}_{}_{}_{}", current_chrom, f_start, f_end, current_type)
                    };
                    features.push(FeatureRecord {
                        chrom: current_chrom.clone(),
                        feature_start: f_start,
                        feature_end: f_end,
                        strand,
                        feature_id: fid,
                        feature_type: std::mem::take(&mut current_type),
                    });
                } else {
                    current_type.clear();
                    current_id.clear();
                }
            }
            in_features = false;
        } else if in_features {
            if line.len() >= 21 && !line[5..21].trim().is_empty() && !line[5..21].contains('/') {
                if !current_type.is_empty() {
                    if let Some((f_start, f_end, strand)) = parse_genbank_location(&current_loc) {
                        let fid = if !current_id.is_empty() {
                            std::mem::take(&mut current_id)
                        } else {
                            format!("{}_{}_{}_{}", current_chrom, f_start, f_end, current_type)
                        };
                        features.push(FeatureRecord {
                            chrom: current_chrom.clone(),
                            feature_start: f_start,
                            feature_end: f_end,
                            strand,
                            feature_id: fid,
                            feature_type: std::mem::take(&mut current_type),
                        });
                    }
                    current_type.clear();
                    current_id.clear();
                    current_loc.clear();
                }

                current_type = line[5..21].trim().to_string();
                current_loc = line[21..].trim().to_string();
            } else if !current_type.is_empty() {
                let trimmed = line.trim();
                if trimmed.starts_with('/') {
                    if let Some((key, val)) = trimmed[1..].split_once('=') {
                        let k = key.trim();
                        let clean_v = val.trim().trim_matches('"').trim();
                        let is_locus_tag = k.eq_ignore_ascii_case("locus_tag");
                        let is_gene_id = k.eq_ignore_ascii_case("gene")
                            || k.eq_ignore_ascii_case("protein_id")
                            || k.eq_ignore_ascii_case("db_xref")
                            || k.eq_ignore_ascii_case("ID");

                        if is_locus_tag || (is_gene_id && current_id.is_empty()) {
                            current_id = clean_v.to_string();
                        }
                    }
                } else if !current_loc.is_empty() && !trimmed.starts_with('/') {
                    current_loc.push_str(trimmed);
                }
            }
        }
    }

    Ok(features)
}

fn parse_genbank_location(loc_str: &str) -> Option<(u32, u32, bool)> {
    let strand = !loc_str.contains("complement");

    let mut numbers = Vec::new();
    let mut current_num = String::new();

    for b in loc_str.bytes() {
        if b.is_ascii_digit() {
            current_num.push(b as char);
        } else if !current_num.is_empty() {
            if let Ok(n) = current_num.parse::<u32>() {
                numbers.push(n);
            }
            current_num.clear();
        }
    }
    if !current_num.is_empty() {
        if let Ok(n) = current_num.parse::<u32>() {
            numbers.push(n);
        }
    }

    if numbers.is_empty() {
        return None;
    }

    let min_pos = *numbers.iter().min()?;
    let max_pos = *numbers.iter().max()?;

    if min_pos == 0 {
        return None;
    }

    let feature_start = min_pos - 1;
    let feature_end = max_pos;

    Some((feature_start, feature_end, strand))
}

/// Read FeatureRecord list transparently from GFF/GTF or GenBank file
pub fn read_feature_records(path: &Path) -> Result<Vec<FeatureRecord>> {
    let mut reader = super::open_compressed_reader(path)?;

    let mut header_buf = [0u8; 1024];
    let n = reader.read(&mut header_buf)?;

    let peek_str = String::from_utf8_lossy(&header_buf[..n]);

    let combined_reader = std::io::Cursor::new(header_buf[..n].to_vec()).chain(reader);

    if peek_str.contains("LOCUS") || peek_str.contains("FEATURES") || peek_str.contains("ORIGIN") {
        parse_genbank_features(combined_reader)
    } else {
        parse_gff_gtf_features(combined_reader)
    }
}

/// Build Polars features DataFrame from FeatureRecord list
pub fn build_features_dataframe(features: &[FeatureRecord]) -> Result<DataFrame> {
    let mut chrom_vec = Vec::with_capacity(features.len());
    let mut start_vec = Vec::with_capacity(features.len());
    let mut end_vec = Vec::with_capacity(features.len());
    let mut strand_vec = Vec::with_capacity(features.len());
    let mut id_vec = Vec::with_capacity(features.len());
    let mut type_vec = Vec::with_capacity(features.len());

    for f in features {
        chrom_vec.push(f.chrom.as_str());
        start_vec.push(f.feature_start);
        end_vec.push(f.feature_end);
        strand_vec.push(f.strand);
        id_vec.push(f.feature_id.as_str());
        type_vec.push(f.feature_type.as_str());
    }

    let chrom_series = Series::new("chrom".into(), chrom_vec)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
    let type_series = Series::new("feature_type".into(), type_vec)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;

    let df = DataFrame::new(vec![
        chrom_series.into(),
        Series::new("feature_start".into(), start_vec).into(),
        Series::new("feature_end".into(), end_vec).into(),
        Series::new("strand".into(), strand_vec).into(),
        Series::new("feature_id".into(), id_vec).into(),
        type_series.into(),
    ])?;

    Ok(df)
}
