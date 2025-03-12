use crate::constants;
use base64::{Engine as _, prelude::BASE64_URL_SAFE};
use itertools::Itertools;
use sha2::{Digest, Sha512};

#[inline(always)]
pub fn sha512t24u_digest_default(attr: &[u8]) -> String {
    sha512t24u_digest(attr, constants::DEFAULT_DIGEST_BYTES)
}

#[inline(always)]
pub fn sha512t24u_digest(attr: &[u8], offset: usize) -> String {
    let digest = Sha512::digest(attr);
    BASE64_URL_SAFE.encode(&digest[..offset])
}

#[inline(always)]
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
