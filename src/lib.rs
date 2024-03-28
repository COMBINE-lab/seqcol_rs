#[derive(Debug)]
struct SeqCol {
    lengths: Vec<usize>,
    names: Vec<String>,
    sequences: Vec<String>
}

impl SeqCol {
    pub fn from_sam_header<'a, I>(it: I) -> Self
        where 
        I: IntoIterator<Item=(&'a [u8], usize)>
    {
        let mut names = Vec::new();
        let mut lengths= Vec::new();
        for (n, l) in it {
            names.push(std::str::from_utf8(n).unwrap().to_owned());
            lengths.push(l);
        }

        Self {
            lengths,
            names,
            sequences: vec![]
        }
    }
}


pub fn add(left: usize, right: usize) -> usize {
    left + right
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
            header.reference_sequences().iter().map(|(k, v)| (k.as_slice(), v.length().into() )) 
        );

        println!("s = {s:?}");
    }


    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
