use base64::{prelude::BASE64_URL_SAFE, Engine as _};
use itertools::Itertools;
use sha2::{Digest, Sha512};

pub fn sha512t24u_digest(attr: &[u8], offset: usize) -> String {
    let digest = Sha512::digest(attr);
    BASE64_URL_SAFE.encode(&digest[..offset])
}

pub fn canonical_rep(inp: &serde_json::Value) -> anyhow::Result<String> {
    if let Some(om) = inp.as_object() {
        let mut output_obj = vec![];
        for k in om.keys().sorted() {
            output_obj.push((k, om[k].clone()))
        }
        let output_val: serde_json::Value = output_obj.into_iter().collect();
        Ok(serde_json::to_string(&output_val)?)
    } else {
        Ok(serde_json::to_string(inp)?)
    }
}
