use anyhow::Context;
use needletail::Sequence;
use serde_json::json;
use sha2::Digest as ShaDigest;
use sha2::Sha256;
use std::path::Path;

pub mod constants;
pub mod utils;

/// holds information relevant to seqcol
/// signatures (and sha256 signatures)
#[derive(Debug)]
pub struct SeqCol {
    lengths: Vec<usize>,
    names: Vec<String>,
    sequences: Option<Vec<String>>,
    sha256_names: Option<String>,
    sha256_seqs: Option<String>,
}

/// trait that converts a DigestResult to
/// a `json` object.
pub trait DigestToJson {
    fn to_json(&self) -> serde_json::Value;
}

/// just a top-level seqcol digest
#[derive(Debug)]
pub struct Level0Digest {
    pub digest: String,
}

impl DigestToJson for Level0Digest {
    fn to_json(&self) -> serde_json::Value {
        let mut repr = serde_json::map::Map::new();
        let lvl: serde_json::Value = 0_u32.into();
        repr.insert("level".to_owned(), lvl);
        repr.insert(
            "seqcol_digest".to_owned(),
            serde_json::Value::String(self.digest.clone()),
        );
        serde_json::json!(repr)
    }
}

/// a level 1 seqcol digest
#[derive(Debug)]
pub struct Level1Digest {
    pub lengths: String,
    pub names: String,
    pub sequences: Option<String>,
    pub sorted_name_length_pairs: Option<String>,
}

impl DigestToJson for Level1Digest {
    fn to_json(&self) -> serde_json::Value {
        let mut repr = serde_json::map::Map::new();
        let lvl: serde_json::Value = 1_u32.into();
        repr.insert("level".to_owned(), lvl);
        repr.insert(
            "lengths".to_owned(),
            serde_json::Value::String(self.lengths.clone()),
        );
        repr.insert(
            "names".to_owned(),
            serde_json::Value::String(self.names.clone()),
        );
        if let Some(seq) = &self.sequences {
            repr.insert(
                "sequences".to_owned(),
                serde_json::Value::String(seq.clone()),
            );
        }
        if let Some(snlp) = &self.sorted_name_length_pairs {
            repr.insert(
                "sorted_name_length_pairs".to_owned(),
                serde_json::Value::String(snlp.clone()),
            );
        }
        serde_json::json!(repr)
    }
}

/// a level 2 seqcol digest
#[derive(Debug)]
pub struct Level2Digest {
    pub lengths: Vec<usize>,
    pub names: Vec<String>,
    pub sequences: Option<Vec<String>>,
    pub sorted_name_length_pairs: Option<String>,
}

impl DigestToJson for Level2Digest {
    fn to_json(&self) -> serde_json::Value {
        if self.sequences.is_some() {
            serde_json::json!({
                "level" : 2,
                "lengths" : self.lengths,
                "names" : self.names,
                "sequences" : self.sequences.clone().expect("non-empty")
            })
        } else {
            serde_json::json!({
                "level" : 2,
                "lengths" : self.lengths,
                "names" : self.names,
            })
        }
    }
}

#[derive(Debug)]
pub enum DigestLevelResult {
    Level0(Level0Digest),
    Level1(Level1Digest),
    Level2(Level2Digest),
}

impl DigestToJson for DigestLevelResult {
    fn to_json(&self) -> serde_json::Value {
        match self {
            Self::Level0(l) => l.to_json(),
            Self::Level1(l) => l.to_json(),
            Self::Level2(l) => l.to_json(),
        }
    }
}

/// Represents the computed result of a digest.
/// This will always have a valid `sq_digest` and
/// may have sha256 digests for the names and or sequences.
#[derive(Debug)]
pub struct DigestResult {
    pub sq_digest: DigestLevelResult,
    pub sha256_names: Option<String>,
    pub sha256_seqs: Option<String>,
}

impl DigestResult {
    /// Produce a [serde_json::Value] JSON representation of the
    /// current [DigestResult].
    pub fn to_json(&self) -> serde_json::Value {
        let mut repr = serde_json::map::Map::new();
        repr.insert("seqcol_digest".to_owned(), self.sq_digest.to_json());
        if self.sha256_names.is_some() || self.sha256_seqs.is_some() {
            repr.insert(
                "sha256_digests".to_owned(),
                json!({
                    "sha256_names" : self.sha256_names,
                    "sha256_seqs" : self.sha256_seqs }),
            );
        } else {
            repr.insert("sha256_digests".to_owned(), serde_json::Value::Null);
        }
        serde_json::json!(repr)
    }
}

#[allow(dead_code)]
pub enum DigestKeyType {
    Required(String),
    Optional(String),
}

#[derive(Debug, Clone, Copy)]
pub enum DigestLevel {
    Level0,
    Level1,
    Level2,
}

#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
/// The configuration describing how a digest should
/// be computed.
pub struct DigestConfig {
    pub level: DigestLevel,
    /// Compute the digest of the
    /// sorted list of (sequence name, length) pair digests
    pub with_seqname_pairs: bool,
}

/// The default configuration uses only the required fields
impl Default for DigestConfig {
    fn default() -> Self {
        Self {
            level: DigestLevel::Level1,
            with_seqname_pairs: false,
        }
    }
}

impl DigestConfig {
    pub fn with_level(level: DigestLevel) -> Self {
        Self {
            level,
            with_seqname_pairs: false,
        }
    }

    pub fn with_level_and_seqname_pairs(level: DigestLevel) -> Self {
        Self {
            level,
            with_seqname_pairs: true,
        }
    }
}

pub struct DigestFunction {
    digest: fn(&[u8]) -> String,
}

impl DigestFunction {
    #[inline(always)]
    pub fn compute(&self, x: &[u8]) -> String {
        (self.digest)(x)
    }
}

impl Default for DigestFunction {
    fn default() -> Self {
        Self {
            digest: utils::sha512t24u_digest_default,
        }
    }
}

#[allow(dead_code)]
impl SeqCol {
    /// Takes an iterator `it` over pairs of sequence names and lengths, and populates
    /// this [SeqCol] object from the contents of the iterator.
    pub fn from_sam_header<'a, I>(it: I) -> Self
    where
        I: IntoIterator<Item = (&'a [u8], usize)>,
    {
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        for (n, l) in it {
            names.push(std::str::from_utf8(n).unwrap().to_owned());
            lengths.push(l);
        }

        Self {
            lengths,
            names,
            sequences: None,
            sha256_names: None,
            sha256_seqs: None,
        }
    }

    /// Tries to construct a [SeqCol] object from an input seq_col structure
    /// represented as a [serde_json::Value] (i.e. deserialized from a JSON file).
    /// This returns an [Ok]`(`[SeqCol]``)` if successful, or an error if the
    /// object couldn't be constructed (e.g. required fields were missing, etc.).
    pub fn try_from_seqcol(sc: &serde_json::Value) -> anyhow::Result<Self> {
        // The  input seqcol object must be of type `Object`
        if let Some(seqcol) = sc.as_object() {
            // we *must* have a lengths field
            let lengths = seqcol
                .get("lengths")
                .ok_or(anyhow::anyhow!("must contain a \"lengths\" field"))?;
            let lengths = lengths
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect();

            // we *must* have a names field
            let names = seqcol
                .get("names")
                .ok_or(anyhow::anyhow!("must contain a \"names\" field"))?;
            let names = names
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_str().unwrap().to_owned())
                .collect();

            // we *may* have a sequences field
            let sequences = seqcol.get("sequences").map(|seqs| {
                seqs.as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_str().unwrap().to_owned())
                    .collect()
            });

            Ok(Self {
                lengths,
                names,
                sequences,
                sha256_names: None,
                sha256_seqs: None,
            })
        } else {
            anyhow::bail!("seqcol object must be a valid JSON Object");
        }
    }

    /// Try to construct a [SeqCol] object from an input FASTA file (represented
    /// by it's [AsRef<Path>] `fp`).
    /// Returns [Ok]`(`[SeqCol]`)` on success or otherwise an error describing
    /// why the [SeqCol] object could not be constructed.
    pub fn try_from_fasta_file<P: AsRef<Path>>(fp: P) -> anyhow::Result<Self> {
        let mut reader = needletail::parse_fastx_file(&fp).with_context(|| {
            format!("cannot parse FASTA records from {}", &fp.as_ref().display())
        })?;

        let digest_function = DigestFunction::default();
        let mut name_sha_digest = Sha256::new();
        let mut seq_sha_digest = Sha256::new();

        let mut names = vec![];
        let mut lengths = vec![];
        let mut seqs = vec![];

        while let Some(record) = reader.next() {
            let seqrec = record?;
            let seq_bytes = seqrec.normalize(false);
            seq_sha_digest.update(seq_bytes.as_ref());
            let h = digest_function.compute(seq_bytes.as_ref());
            seqs.push(format!("SQ.{h}"));

            // take the record name up to the first whitespace
            if let Some(seq_name) = std::str::from_utf8(seqrec.id())?
                .split_whitespace()
                .next()
                .map(str::to_owned)
            {
                name_sha_digest.update(seq_name.as_bytes());
                names.push(seq_name);
            } else {
                anyhow::bail!("cannot process data with empty sequence names!")
            };

            lengths.push(seqrec.num_bases());
        }

        let sha256_names = hex::encode(name_sha_digest.finalize());
        let sha256_seqs = hex::encode(seq_sha_digest.finalize());

        Ok(Self {
            lengths,
            names,
            sequences: Some(seqs),
            sha256_names: Some(sha256_names),
            sha256_seqs: Some(sha256_seqs),
        })
    }

    /// Takes an iterator `it` over pairs of sequence names and seqeuences, and populates
    /// this [SeqCol] object from the contents of the iterator.
    pub fn try_from_name_seq_iter<I>(it: I) -> anyhow::Result<Self>
    where
        I: IntoIterator<Item = (String, String)>,
    {
        let digest_function = DigestFunction::default();
        let mut name_sha_digest = Sha256::new();
        let mut seq_sha_digest = Sha256::new();

        let mut names = vec![];
        let mut lengths = vec![];
        let mut seqs = vec![];

        // iterate over the name sequence pairs and add them
        // to the seqcol digest (and update the sha digests)
        for (name, seq) in it {
            seq_sha_digest.update(&seq);
            let h = digest_function.compute(seq.as_ref());
            seqs.push(format!("SQ.{h}"));

            name_sha_digest.update(name.as_bytes());
            names.push(name);

            lengths.push(seq.len());
        }

        let sha256_names = hex::encode(name_sha_digest.finalize());
        let sha256_seqs = hex::encode(seq_sha_digest.finalize());

        Ok(Self {
            lengths,
            names,
            sequences: Some(seqs),
            sha256_names: Some(sha256_names),
            sha256_seqs: Some(sha256_seqs),
        })
    }

    /// Computes and returns the [SeqCol] digest of the current [SeqCol] object.
    /// The [DigestConfig] parameter `c` controls what fields are used to compute the digest.
    ///
    /// Returns [Ok]`(`[String]`)` on success, representing the computed digest
    /// or otherwise an error describing why the digest could not be computed.
    pub fn digest(&self, c: DigestConfig) -> anyhow::Result<DigestResult> {
        let digest_function = DigestFunction::default();
        match c.level {
            DigestLevel::Level0 => {
                let sq_json = self.seqcol_obj(c)?;
                let mut digest_json = json!({});
                for (k, v) in sq_json.as_object().unwrap().iter() {
                    let v2 = utils::canonical_rep(v)?;
                    let h2 = digest_function.compute(v2.as_bytes());
                    digest_json[k] = serde_json::Value::String(h2);
                }
                let digest_str = utils::canonical_rep(&digest_json)?;
                Ok(DigestResult {
                    sq_digest: DigestLevelResult::Level0(Level0Digest {
                        digest: digest_function.compute(digest_str.as_bytes()),
                    }),
                    sha256_names: self.sha256_names.clone(),
                    sha256_seqs: self.sha256_seqs.clone(),
                })
            }
            DigestLevel::Level1 => {
                let sq_json = self.seqcol_obj(c)?;
                let mut names = String::new();
                let mut lengths = String::new();
                let mut sequences = None;
                let mut sorted_name_length_pairs = None;
                for (k, v) in sq_json.as_object().unwrap().iter() {
                    let v2 = utils::canonical_rep(v)?;
                    let h2 = digest_function.compute(v2.as_bytes());
                    match k.as_str() {
                        "names" => {
                            names = h2;
                        }
                        "lengths" => {
                            lengths = h2;
                        }
                        "sequences" => {
                            sequences = Some(h2);
                        }
                        "sorted_name_length_pairs" => {
                            sorted_name_length_pairs = Some(h2);
                        }
                        _ => {}
                    }
                }
                Ok(DigestResult {
                    sq_digest: DigestLevelResult::Level1(Level1Digest {
                        names,
                        lengths,
                        sequences,
                        sorted_name_length_pairs,
                    }),
                    sha256_names: self.sha256_names.clone(),
                    sha256_seqs: self.sha256_seqs.clone(),
                })
            }
            DigestLevel::Level2 => {
                let snlp_digest_val = if c.with_seqname_pairs {
                    let pairs = serde_json::Value::Array(self.get_sorted_name_len_pair_digests()?);
                    let v = utils::canonical_rep(&pairs)?;
                    Some(digest_function.compute(v.as_bytes()))
                } else {
                    None
                };
                Ok(DigestResult {
                    sq_digest: DigestLevelResult::Level2(Level2Digest {
                        names: self.names.clone(),
                        lengths: self.lengths.clone(),
                        sequences: self.sequences.clone(),
                        sorted_name_length_pairs: snlp_digest_val,
                    }),
                    sha256_names: self.sha256_names.clone(),
                    sha256_seqs: self.sha256_seqs.clone(),
                })
            }
        }
    }

    fn get_sorted_name_len_pair_digests(&self) -> anyhow::Result<Vec<serde_json::Value>> {
        let mut snlp_digests = vec![];
        snlp_digests.reserve_exact(self.names.len());

        let digest_function = DigestFunction::default();
        for (l, n) in self.lengths.iter().zip(self.names.iter()) {
            let cr = utils::canonical_rep(&json!({"length" : l, "name" : n}))?;
            snlp_digests.push(digest_function.compute(cr.as_bytes()));
        }
        snlp_digests.sort_unstable();
        Ok(snlp_digests
            .into_iter()
            .map(serde_json::Value::String)
            .collect())
    }

    pub fn seqcol_obj(&self, c: DigestConfig) -> anyhow::Result<serde_json::Value> {
        let mut sq_json = json!({
            "lengths" : self.lengths,
            "names" : self.names,
        });

        if let Some(v) = &self.sequences {
            sq_json["sequences"] = serde_json::Value::Array(
                v.iter()
                    .map(|x| serde_json::Value::String(x.to_string()))
                    .collect(),
            );
        };

        if c.with_seqname_pairs {
            sq_json["sorted_name_length_pairs"] =
                serde_json::Value::Array(self.get_sorted_name_len_pair_digests()?);
        }
        Ok(sq_json)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn from_sam_header_works() {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tndls
";
        let header: noodles_sam::Header = s.parse().unwrap();

        let s = SeqCol::from_sam_header(
            header
                .reference_sequences()
                .iter()
                .map(|(k, v)| (k.as_slice(), v.length().into())),
        );
        let r = s
            .digest(DigestConfig {
                level: DigestLevel::Level0,
                with_seqname_pairs: false,
            })
            .unwrap();
        match r.sq_digest {
            DigestLevelResult::Level0(l0r) => {
                assert_eq!(l0r.digest, "2HqWKZw8F4VY7q9sfYRM-JJ_RaMXv1eK")
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn from_seqcol_object() {
        let file = File::open("test_data/seqcol_obj.json").expect("can't open input seqcol file");
        let reader = BufReader::new(file);
        let sc = serde_json::from_reader(reader).unwrap();
        let s = SeqCol::try_from_seqcol(&sc).unwrap();
        let r = s
            .digest(DigestConfig {
                level: DigestLevel::Level0,
                with_seqname_pairs: false,
            })
            .unwrap();
        match r.sq_digest {
            DigestLevelResult::Level0(l0r) => {
                assert_eq!(l0r.digest, "2HqWKZw8F4VY7q9sfYRM-JJ_RaMXv1eK")
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn from_fasta_file_works_with_default() {
        let s = SeqCol::try_from_fasta_file(Path::new("test_data/simple.fa")).unwrap();
        let r = s
            .digest(DigestConfig {
                level: DigestLevel::Level0,
                with_seqname_pairs: false,
            })
            .unwrap();
        match r.sq_digest {
            DigestLevelResult::Level0(l0r) => {
                assert_eq!(l0r.digest, "E0cJxnAB5lrWXGP_JoWRNWKEDfdPUDUR");
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn from_fasta_file_works_with_default2() {
        let s = SeqCol::try_from_fasta_file(Path::new("test_data/simple2.fa")).unwrap();
        let r = s
            .digest(DigestConfig {
                level: DigestLevel::Level0,
                with_seqname_pairs: false,
            })
            .unwrap();
        match r.sq_digest {
            DigestLevelResult::Level0(l0r) => {
                assert_eq!(l0r.digest, "zsWu-iN7EJt5-8_inZOugaAg3eT-unK3");
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn from_fasta_file_works_with_seqname_length_pairs() {
        let s = SeqCol::try_from_fasta_file(Path::new("test_data/simple.fa")).unwrap();
        let r = s
            .digest(DigestConfig {
                level: DigestLevel::Level0,
                with_seqname_pairs: true,
            })
            .unwrap();
        match r.sq_digest {
            DigestLevelResult::Level0(l0r) => {
                assert_eq!(l0r.digest, "bXpsYPctlKYGMvDGwmoHTUuS7ryH5miY");
            }
            _ => unreachable!(),
        }
        assert_eq!(
            r.sha256_names,
            Some("d3a44a65dac11f07a2aeea6f505e8134dad9e2f9af9c162b83f6df0ce6adbbe5".to_owned())
        );
        assert_eq!(
            r.sha256_seqs,
            Some("8a4b348914002bcc60c8bb2a938d0eef2bb519bd575ec12536d39208e4d3ea4a".to_owned())
        );
    }
}
