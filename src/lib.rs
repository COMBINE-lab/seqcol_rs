use anyhow::Context;
use needletail::Sequence;
use serde_json::json;
use sha2::Digest as ShaDigest;
use sha2::Sha256;
use std::path::Path;

pub mod constants;
pub mod utils;

pub(crate) const INHERENT_ATTRIBUTES: [&str; 2] = ["names", "sequences"];

/// holds information relevant to seqcol
/// signatures (and sha256 signatures)
#[derive(Debug)]
pub struct SeqCol {
    attributes: Vec<SeqColAttribute>,
    sha256_names: Option<String>,
    sha256_seqs: Option<String>,
}

impl SeqCol {
    pub fn names(&self) -> anyhow::Result<Vec<String>> {
        let mut v = None;
        for x in &self.attributes {
            if let SeqColAttribute::Names(x) = x {
                v = Some(x.clone());
            };
        }
        v.ok_or(anyhow::anyhow!(
            "The attribute 'names' was not present in this seqcol object."
        ))
    }

    pub fn lengths(&self) -> anyhow::Result<Vec<usize>> {
        let mut v = None;
        for x in &self.attributes {
            if let SeqColAttribute::Lengths(x) = x {
                v = Some(x.clone());
            }
        }
        v.ok_or(anyhow::anyhow!(
            "The attribute 'length' was not present in this seqcol object."
        ))
    }

    pub fn sequences(&self) -> anyhow::Result<Vec<String>> {
        let mut v = None;
        for x in &self.attributes {
            if let SeqColAttribute::Sequences(x) = x {
                v = Some(x.clone());
            }
        }
        v.ok_or(anyhow::anyhow!(
            "The attribute 'sequences' was not present in this seqcol object."
        ))
    }
}

#[derive(Debug)]
pub enum SeqColAttribute {
    Lengths(Vec<usize>),
    Names(Vec<String>),
    Sequences(Vec<String>),
    SortedSequences(Vec<String>),
    NameLengthPairs(Vec<(String, usize)>),
    SortedNameLengthPairs(Vec<(String, usize)>),
}

impl SeqColAttribute {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Lengths(_) => "lengths",
            Self::Names(_) => "names",
            Self::Sequences(_) => "sequences",
            Self::SortedSequences(_) => "sorted_sequences",
            Self::NameLengthPairs(_) => "name_length_pairs",
            Self::SortedNameLengthPairs(_) => "sorted_name_length_pairs",
        }
    }

    pub fn is_inherent(&self) -> bool {
        matches!(self, Self::Names(_) | Self::Sequences(_))
    }

    pub fn is_required(&self) -> bool {
        matches!(self, Self::Names(_) | Self::Sequences(_) | Self::Lengths(_))
    }

    pub fn is_collated(&self) -> bool {
        matches!(
            self,
            Self::Names(_) | Self::Sequences(_) | Self::Lengths(_) | Self::NameLengthPairs(_)
        )
    }

    fn try_into_level2_repr(&self) -> anyhow::Result<serde_json::Value> {
        match self {
            Self::Lengths(l) => Ok(serde_json::Value::Array(
                l.iter().map(|x| json!(*x as u64)).collect(),
            )),
            Self::Names(n) => Ok(serde_json::Value::Array(
                n.iter().map(|x| json!(x)).collect(),
            )),
            Self::Sequences(s) => Ok(serde_json::Value::Array(
                s.iter().map(|x| json!(x)).collect(),
            )),
            Self::SortedSequences(s) => Ok(serde_json::Value::Array(
                s.iter().map(|x| json!(x)).collect(),
            )),
            Self::NameLengthPairs(nlp) => Ok(serde_json::Value::Array(
                nlp.iter()
                    .map(|(n, l)| json!({ "name" : n, "length" : *l as u64}))
                    .collect(),
            )),
            Self::SortedNameLengthPairs(snlp) => {
                let digest_function = DigestFunction::default();
                let mut digests: Vec<String> = snlp
                    .iter()
                    .map(|(n, l)| {
                        let j = json!({ "name" : n, "length" : *l as u64});
                        let v2 = utils::canonical_rep(&j).expect("should have canoncal repr");
                        let h2 = digest_function.compute(v2.as_bytes());
                        h2
                    })
                    .collect();
                digests.sort_unstable();
                Ok(json!(digests))
            }
        }
    }

    pub fn try_into_level_repr(&self, level: DigestLevel) -> anyhow::Result<serde_json::Value> {
        match level {
            DigestLevel::Level2 => self.try_into_level2_repr(),
            DigestLevel::Level1 => {
                let digest_function = DigestFunction::default();
                let l2 = self.try_into_level2_repr()?;
                let v2 = utils::canonical_rep(&l2)?;
                Ok(serde_json::Value::String(
                    digest_function.compute(v2.as_bytes()),
                ))
            }
            DigestLevel::Level0 => {
                anyhow::bail!("It does not make sense to try to convert an individual seqcol attribute into a level 0 representation")
            }
        }
    }
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
            "digest".to_owned(),
            serde_json::Value::String(self.digest.clone()),
        );
        serde_json::json!(repr)
    }
}

/// a level 1 seqcol digest
#[derive(Debug)]
pub struct Level1Digest {
    digests: serde_json::Value,
}

impl DigestToJson for Level1Digest {
    fn to_json(&self) -> serde_json::Value {
        let mut repr = serde_json::map::Map::new();
        let lvl: serde_json::Value = 1_u32.into();
        repr.insert("level".to_owned(), lvl);
        for (k, v) in self
            .digests
            .as_object()
            .expect("should be an object")
            .clone()
        {
            repr.insert(k, v);
        }

        serde_json::json!(repr)
    }
}

/// a level 2 seqcol digest
#[derive(Debug)]
pub struct Level2Digest {
    pub digests: serde_json::Value,
}

impl DigestToJson for Level2Digest {
    fn to_json(&self) -> serde_json::Value {
        let mut repr = serde_json::map::Map::new();
        let lvl: serde_json::Value = 2_u32.into();
        repr.insert("level".to_owned(), lvl);
        for (k, v) in self
            .digests
            .as_object()
            .expect("should be an object")
            .clone()
        {
            repr.insert(k, v);
        }
        serde_json::json!(repr)
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

#[derive(Debug, Clone, Copy)]
pub enum DigestLevel {
    Level0,
    Level1,
    Level2,
}

#[derive(Debug, Clone, Copy)]
pub enum KnownAttr {
    NameLengthPairs,
    SortedNameLengthPairs,
    SortedSequences,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
/// The configuration describing how a digest should
/// be computed.
pub struct DigestConfig {
    pub level: DigestLevel,
    /// additional attributes to include
    pub additional_attr: Vec<KnownAttr>,
}

/// The default configuration uses only the required fields
impl Default for DigestConfig {
    fn default() -> Self {
        Self {
            level: DigestLevel::Level1,
            additional_attr: vec![],
        }
    }
}

impl DigestConfig {
    pub fn with_level(level: DigestLevel) -> Self {
        Self {
            level,
            additional_attr: vec![],
        }
    }

    pub fn with_level_and_additional_attrs(level: DigestLevel, attr: &[KnownAttr]) -> Self {
        Self {
            level,
            additional_attr: attr.to_vec(),
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

        let attributes = vec![
            SeqColAttribute::Names(names),
            SeqColAttribute::Lengths(lengths),
        ];

        Self {
            attributes,
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

            let mut attributes = vec![
                SeqColAttribute::Names(names),
                SeqColAttribute::Lengths(lengths),
            ];
            if let Some(seqs) = sequences {
                attributes.push(SeqColAttribute::Sequences(seqs));
            }
            Ok(Self {
                attributes,
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

        let attributes = vec![
            SeqColAttribute::Names(names),
            SeqColAttribute::Lengths(lengths),
            SeqColAttribute::Sequences(seqs),
        ];
        Ok(Self {
            attributes,
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

        let attributes = vec![
            SeqColAttribute::Names(names),
            SeqColAttribute::Lengths(lengths),
            SeqColAttribute::Sequences(seqs),
        ];
        Ok(Self {
            attributes,
            sha256_names: Some(sha256_names),
            sha256_seqs: Some(sha256_seqs),
        })
    }

    fn get_derived_attr(&self, att: KnownAttr) -> anyhow::Result<SeqColAttribute> {
        match att {
            KnownAttr::NameLengthPairs => {
                let names = self.names()?;
                let lens = self.lengths()?;
                assert_eq!(names.len(), lens.len());
                Ok(SeqColAttribute::NameLengthPairs(
                    names.into_iter().zip(lens).collect(),
                ))
            }
            KnownAttr::SortedNameLengthPairs => {
                let names = self.names()?;
                let lens = self.lengths()?;
                assert_eq!(names.len(), lens.len());
                let mut nlp: Vec<(String, usize)> = names.into_iter().zip(lens).collect();
                nlp.sort_unstable();
                Ok(SeqColAttribute::SortedNameLengthPairs(nlp))
            }
            KnownAttr::SortedSequences => {
                let mut sequences = self.sequences()?;
                sequences.sort_unstable();
                Ok(SeqColAttribute::SortedSequences(sequences))
            }
        }
    }

    /// Computes and returns the [SeqCol] digest of the current [SeqCol] object.
    /// The [DigestConfig] parameter `c` controls what fields are used to compute the digest.
    ///
    /// Returns [Ok]`(`[String]`)` on success, representing the computed digest
    /// or otherwise an error describing why the digest could not be computed.
    pub fn digest(&mut self, c: DigestConfig) -> anyhow::Result<DigestResult> {
        let digest_function = DigestFunction::default();
        match c.level {
            DigestLevel::Level0 => {
                let mut inherent_set: Vec<&str> = vec![];
                let mut digest_json = json!({});
                for attr in &self.attributes {
                    if attr.is_inherent() {
                        digest_json[attr.name()] = attr.try_into_level_repr(DigestLevel::Level1)?;
                        inherent_set.push(attr.name());
                    } else if !attr.is_required() {
                        eprintln!("Note: The level 0 digest depends only on inherent attributes, but {} is neither inherent nor required. \
                                   Consider not requesting this attribute at all if you only want a level 0 digest", attr.name());
                    }
                }
                inherent_set.sort_unstable();
                if inherent_set.len() != INHERENT_ATTRIBUTES.len() {
                    anyhow::bail!("level 0 digest requires all inherent attributes {:?}, but this seqcol object had only {:?}", INHERENT_ATTRIBUTES, inherent_set);
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
                let mut additional_attr = Vec::with_capacity(c.additional_attr.len());
                for requested_attr in c.additional_attr {
                    additional_attr.push(self.get_derived_attr(requested_attr)?);
                }
                self.attributes.extend(additional_attr);

                let mut res = json!({});
                for attr in &self.attributes {
                    res[attr.name()] = attr.try_into_level_repr(DigestLevel::Level1)?;
                }
                Ok(DigestResult {
                    sq_digest: DigestLevelResult::Level1(Level1Digest { digests: res }),
                    sha256_names: self.sha256_names.clone(),
                    sha256_seqs: self.sha256_seqs.clone(),
                })
            }
            DigestLevel::Level2 => {
                let mut res = json!({});

                let mut additional_attr = Vec::with_capacity(c.additional_attr.len());
                for requested_attr in c.additional_attr {
                    additional_attr.push(self.get_derived_attr(requested_attr)?);
                }
                self.attributes.extend(additional_attr);

                for attr in &self.attributes {
                    res[attr.name()] = attr.try_into_level_repr(DigestLevel::Level2)?;
                }
                Ok(DigestResult {
                    sq_digest: DigestLevelResult::Level2(Level2Digest { digests: res }),
                    sha256_names: self.sha256_names.clone(),
                    sha256_seqs: self.sha256_seqs.clone(),
                })
            }
        }
    }

    pub fn seqcol_obj(&mut self, c: DigestConfig) -> anyhow::Result<serde_json::Value> {
        let mut sq_json = json!({});

        let mut additional_attr = Vec::with_capacity(c.additional_attr.len());
        for requested_attr in c.additional_attr {
            additional_attr.push(self.get_derived_attr(requested_attr)?);
        }
        self.attributes.extend(additional_attr);

        for attr in &self.attributes {
            sq_json[attr.name()] = attr.try_into_level_repr(DigestLevel::Level2)?;
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
