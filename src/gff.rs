use anyhow::Result;
use noodles::gff;
use rust_lapper::Interval;
use std::fs::File;
use std::io::BufReader;

type Iv = Interval<usize, ()>;

#[derive(Debug)]
pub struct Gene {
    pub name: String,
    pub exons: Vec<Iv>,
    pub _introns: Vec<Iv>,
    pub tid: String,
    pub _start: usize,
    pub _end: usize,
}

pub(crate) fn parse_gff3(gff3_file: &str) -> Result<(Vec<Gene>, usize)> {
    let mut reader = File::open(gff3_file)
        .map(BufReader::new)
        .map(gff::io::Reader::new)?;

    let mut gene_result = Vec::new();
    let mut gene_count = 0;

    for result in reader.records() {
        let record = result?;

        if record.ty() != "gene" && record.ty() != "exon" {
            continue;
        } else {
            // generate a Gene when we see a gene record
            // push to gene_result when we see a new gene record
            if record.ty() == "gene" {
                let gene = Gene {
                    name: record.attributes().get("Name").unwrap().to_string(),
                    exons: Vec::new(),
                    _introns: Vec::new(),
                    tid: record.reference_sequence_name().to_string(),
                    _start: record.start().get(),
                    _end: record.end().get(),
                };
                gene_result.push(gene);
                gene_count += 1;
            } else {
                // add exon to the last gene in gene_result
                let exon_iv = Iv {
                    start: record.start().get(),
                    stop: record.end().get(),
                    val: (),
                };
                gene_result.last_mut().unwrap().exons.push(exon_iv);
            }
        }
    }

    Ok((gene_result, gene_count))
}
