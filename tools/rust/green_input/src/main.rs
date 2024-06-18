fn main() {
    println!("green_input in RUST");
    for token in xmlparser::Tokenizer::from("<tagname name='value'/>") {
        println!("{:?}", token);
    }
}
