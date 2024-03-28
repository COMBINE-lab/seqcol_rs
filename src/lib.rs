use serde_json::json;

pub mod utils;

#[derive(Debug)]
struct SeqCol {
    lengths: Vec<usize>,
    names: Vec<String>,
    sequences: Option<Vec<String>>,
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

    pub fn digest(&self) -> anyhow::Result<String> {
        let empty = Vec::<String>::new();
        let sq_json = json!({
            "lengths" : self.lengths,
            "names" : self.names,
            "sequences" : match &self.sequences {
                Some(v) => v,
                None => &empty
            }
        });

        let mut digest_json = json!({});
        for (k, v) in sq_json.as_object().unwrap().iter() {
            let v2 = utils::canonical_rep(v)?;
            let h2 = utils::sha512t24u_digest(v2.as_bytes(), 24);
            digest_json[k] = serde_json::Value::String(h2);
        }
        let digest_str = utils::canonical_rep(&digest_json)?;
        Ok(utils::sha512t24u_digest(digest_str.as_bytes(), 24))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam;

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
        let r = s.digest().unwrap();

        assert_eq!(r, "6_Sn0CtEZ-LIJDPyhIwYQFBEFnAxDE2j");
    }
}
