use anyhow::Context;
use needletail::Sequence;
use serde_json::json;
use std::path::Path;

pub mod constants;
pub mod utils;

#[derive(Debug)]
pub struct SeqCol {
    lengths: Vec<usize>,
    names: Vec<String>,
    sequences: Option<Vec<String>>,
}

pub enum DigestKeyType {
    Required(String),
    Optional(String),
}

#[allow(dead_code)]
#[derive(Debug)]
/// The configuration describing how a digest should
/// be computed.
pub enum DigestConfig {
    /// Compute the digest using only the required fields
    /// of "lengths", "names" and "sequences" (if present)
    RequiredOnly,
    /// Compute the digest additionally including the
    /// sorted list of (sequence name, length) pair digests
    WithSeqnameLenPairs,
}

/// The default configuration uses only the required fields
impl Default for DigestConfig {
    fn default() -> Self {
        Self::RequiredOnly
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
        }
    }

    pub fn try_from_seqcol(sc: &serde_json::Value) -> anyhow::Result<Self> {
        // The  input seqcol object must be of type `Object`
        if let Some(seqcol) = sc.as_object() {
            // we *must* have a lengths field
            let lengths = seqcol
                .get("lengths")
                .ok_or(anyhow::anyhow!("must contain lengths field"))?;
            let lengths = lengths
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect();

            // we *must* have a names field
            let names = seqcol
                .get("names")
                .ok_or(anyhow::anyhow!("must contain names field"))?;
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
            })
        } else {
            anyhow::bail!("seqcol object must be a valid JSON Object");
        }
    }

    pub fn try_from_fasta_file<P: AsRef<Path>>(fp: P) -> anyhow::Result<Self> {
        let mut reader = needletail::parse_fastx_file(&fp).with_context(|| {
            format!("cannot parse FASTA records from {}", &fp.as_ref().display())
        })?;
        let mut names = vec![];
        let mut lengths = vec![];
        let mut seqs = vec![];
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let h = utils::sha512t24u_digest_default(seqrec.normalize(false).as_ref());
            seqs.push(format!("SQ.{h}"));
            names.push(std::str::from_utf8(seqrec.id())?.to_owned());
            lengths.push(seqrec.num_bases());
        }

        Ok(Self {
            lengths,
            names,
            sequences: Some(seqs),
        })
    }

    pub fn digest(&self, c: DigestConfig) -> anyhow::Result<String> {
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

        let mut snlp_digests = vec![];
        if let DigestConfig::WithSeqnameLenPairs = c {
            snlp_digests.reserve_exact(self.names.len());

            for (l, n) in self.lengths.iter().zip(self.names.iter()) {
                let cr = utils::canonical_rep(&json!({"length" : l, "name" : n}))?;
                snlp_digests.push(utils::sha512t24u_digest_default(cr.as_bytes()));
            }
            snlp_digests.sort_unstable();

            sq_json["sorted_name_length_pairs"] = serde_json::Value::Array(
                snlp_digests
                    .into_iter()
                    .map(serde_json::Value::String)
                    .collect(),
            );
        }

        let mut digest_json = json!({});
        for (k, v) in sq_json.as_object().unwrap().iter() {
            let v2 = utils::canonical_rep(v)?;
            let h2 = utils::sha512t24u_digest_default(v2.as_bytes());
            digest_json[k] = serde_json::Value::String(h2);
        }
        let digest_str = utils::canonical_rep(&digest_json)?;
        Ok(utils::sha512t24u_digest_default(digest_str.as_bytes()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam;
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
        let r = s.digest(DigestConfig::default()).unwrap();

        assert_eq!(r, "2HqWKZw8F4VY7q9sfYRM-JJ_RaMXv1eK");
    }

    #[test]
    fn from_seqcol_object() {
        let file = File::open("test_data/seqcol_obj.json").expect("can't open input seqcol file");
        let reader = BufReader::new(file);
        let sc = serde_json::from_reader(reader).unwrap();
        let s = SeqCol::try_from_seqcol(&sc).unwrap();
        let r = s.digest(DigestConfig::default()).unwrap();
        assert_eq!(r, "2HqWKZw8F4VY7q9sfYRM-JJ_RaMXv1eK");
    }

    #[test]
    fn from_fasta_file_works_with_default() {
        let s = SeqCol::try_from_fasta_file(Path::new("test_data/simple.fa")).unwrap();
        let r = s.digest(DigestConfig::default()).unwrap();
        assert_eq!(r, "E0cJxnAB5lrWXGP_JoWRNWKEDfdPUDUR");
    }

    #[test]
    fn from_fasta_file_works_with_seqname_length_pairs() {
        let s = SeqCol::try_from_fasta_file(Path::new("test_data/simple.fa")).unwrap();
        let r = s.digest(DigestConfig::WithSeqnameLenPairs).unwrap();
        assert_eq!(r, "bXpsYPctlKYGMvDGwmoHTUuS7ryH5miY");
    }
}
