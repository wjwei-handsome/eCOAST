use anyhow::Result;
use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::record::Record;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Cigar, bam::Read};
use std::fmt;
// use smartstring::alias::String;
use std::path::PathBuf;
use std::{convert::TryFrom, rc::Rc};
use std::{default, fmt::Display};

/// A serializable object meant to hold all information about a position.
pub trait Position: Default {
    /// Create a new position with all other values zeroed
    fn new(ref_seq: String, pos: u32) -> Self;
}
/// Hold all information about a range of positions.
#[derive(Debug, Default)]
pub struct BedGraph {
    /// Reference sequence name.
    pub ref_seq: String,
    /// 1-based position in the sequence.
    pub pos: u32,
    /// The point at which this depth ends, non-inclusive
    pub end: u32,
    /// Total depth at this position.
    pub depth: u32,
}

impl Display for BedGraph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}",
            self.ref_seq, self.pos, self.end, self.depth
        )
    }
}

impl Position for BedGraph {
    /// Create a new position for the given ref_seq name.
    fn new(ref_seq: String, pos: u32) -> Self {
        BedGraph {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }
}

/// Bed-like formatting of [`RangePositions`]

/// Anything that implements ReadFilter can apply a filter set to read.
pub trait ReadFilter {
    /// filters a read, true is pass, false if fail
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool;
}

/// A straightforward read filter.
pub struct DefaultReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

impl DefaultReadFilter {
    /// Create an OnlyDepthReadFilter
    pub fn new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self {
        Self {
            include_flags,
            exclude_flags,
            min_mapq,
        }
    }
}

impl ReadFilter for DefaultReadFilter {
    /// Filter reads based SAM flags and mapping quality
    #[inline(always)]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        let flags = read.flags();
        (!flags) & self.include_flags == 0
            && flags & self.exclude_flags == 0
            && read.mapq() >= self.min_mapq
    }
}

// A tweaked impl of IterAlignedBlocks from [here](https://github.com/rust-bio/rust-htslib/blob/9175d3ca186baef4f84a7d7ccb27869b43471e36/src/bam/ext.rs#L51)
// Not that this will also hang onto the bam::Record and supplies the qname for each thing returned.
// At the end of the day this shouldn't be the worst since any given read should not have that many splits in it
struct IterAlignedBlocks {
    pos: i64,
    cigar_index: usize,
    cigar: bam::record::CigarStringView,
    overlap_status: bool,
    record: Rc<bam::Record>,
}
impl IterAlignedBlocks {
    fn new(record: Rc<bam::Record>) -> Self {
        let overlap = false;
        Self {
            pos: record.reference_start(),
            cigar_index: 0,
            cigar: record.cigar(),
            overlap_status: overlap,
            record,
        }
    }
}

impl Iterator for IterAlignedBlocks {
    type Item = (i64, i64, bool, String);
    fn next(&mut self) -> Option<Self::Item> {
        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                    let out_pos = self.pos;
                    self.pos += len as i64;
                    self.cigar_index += 1;
                    return Some((
                        out_pos,
                        out_pos + len as i64,
                        self.overlap_status,
                        String::from(
                            std::str::from_utf8(self.record.qname()).expect("Convert qname"),
                        ),
                    ));
                }
                Cigar::RefSkip(len) => self.pos += len as i64,
                _ => (),
            }
            self.cigar_index += 1;
        }
        None
    }
}

pub(crate) struct OnlyDepthProcessor<F: ReadFilter + Send> {
    /// path to indexed BAM/CRAM
    pub reads: PathBuf,
    /// implementation of [position::ReadFilter] that will be used
    pub read_filter: F,
    /// 0-based or 1-based coordinate output
    pub coord_base: u32,
}

impl<F: ReadFilter + Send> OnlyDepthProcessor<F> {
    /// Create a new OnlyDepthProcessor
    pub fn new(reads: PathBuf, coord_base: u32, read_filter: F) -> Self {
        Self {
            reads,
            coord_base,
            read_filter,
        }
    }

    /// Sum the counts within the region to get the depths at each RangePosition
    #[inline]
    fn sum_counter(
        &self,
        counter: Vec<i32>,
        contig: &str,
        region_start: u32,
        region_stop: u32,
    ) -> Result<Vec<BedGraph>> {
        // Sum the counter and merge same-depth ranges of positions
        let mut sum: i32 = 0;
        let mut results = vec![];
        let mut curr_start = region_start;
        if counter.is_empty() {
            return Ok(results);
        }
        let mut curr_depth = counter[0];
        for (i, count) in counter.iter().enumerate() {
            sum += count;
            // freeze pos and start new one
            if curr_depth != sum {
                let mut pos = BedGraph::new(String::from(contig), curr_start + self.coord_base);
                pos.depth = u32::try_from(curr_depth)?;
                pos.end = region_start + i as u32 + self.coord_base;

                curr_start = region_start + i as u32;
                curr_depth = sum;
                results.push(pos);
            }

            if i == counter.len() - 1 {
                // We've hit the end
                let mut pos = BedGraph::new(String::from(contig), curr_start + self.coord_base);
                pos.depth = u32::try_from(curr_depth)?;
                pos.end = region_stop + self.coord_base;
                results.push(pos);
            }
        }
        Ok(results)
    }

    /// Process a region, taking into account REF_SKIPs and mates
    pub fn process_region(&self, tid: &str, start: u32, stop: u32) -> Result<Vec<BedGraph>> {
        // Create a reader
        let mut reader = bam::IndexedReader::from_path(&self.reads)?;

        // fetch the region of interest
        reader.fetch((tid, start, stop))?;

        let mut counter: Vec<i32> = vec![0; (stop - start) as usize];

        // Walk over each read, counting the starts and ends
        for record in reader
            .rc_records()
            .map(|r| r.unwrap())
            .filter(|read| self.read_filter.filter_read(read, None))
            .flat_map(IterAlignedBlocks::new)
        {
            let rec_start = u32::try_from(record.0)?;
            let rec_stop = u32::try_from(record.1)?;

            // NB: since we are splitting the region, it's possible the region we are looking at
            // may occur before the ROI, or after the ROI
            if rec_start >= stop || rec_stop <= start {
                continue;
            }

            // rectify start / stop with region boundaries
            // increment the start of the region
            let adjusted_start = if rec_start < start {
                0
            } else {
                (rec_start - start) as usize
            };

            let mut dont_count_stop = false; // if this interval extends past the end of the region, don't count an end for it
            let adjusted_stop = if rec_stop >= stop {
                dont_count_stop = true;
                counter.len() - 1
            } else {
                (rec_stop - start) as usize
            };

            // check if this read has a mate that will be seen within this region
            // that this works for both mates in pair
            {
                counter[adjusted_start] += 1;
                if !dont_count_stop {
                    // check if the end of interval extended past region end
                    counter[adjusted_stop] -= 1;
                }
            }
        }

        // Sum the counter and merge same-depth ranges of positions
        // let contig = std::str::from_utf8(header.tid2name(tid)).unwrap();
        self.sum_counter(counter, tid, start, stop)
    }
}
